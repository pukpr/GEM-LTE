using JSON
using Statistics
using Printf

# Ranked wrapper for lte.jl
#
# Usage:
#   julia lte_super.jl <experiments_root> [lambda] [evaluation_root] [max_runs] [--seed-mean <mean_dir>]
#
# Purpose:
#   1. Parse the ranked cross-validation ordering from docs/gem-lte-results.html
#   2. Run lte.jl over all available experiment subdirectories in that order
#   3. Maintain the shared mean manifold through lte.jl's evaluation/mean workflow
#   4. Snapshot the best cumulative mean manifold according to a balanced score
#
# The wrapper does not modify lte.jl behavior; it uses the JSON summaries that
# lte.jl writes into each evaluation directory.

const DEFAULT_RESULTS_HTML = normpath(joinpath(@__DIR__, "..", "gem-lte-results.html"))
const DEFAULT_LTE_SCRIPT = joinpath(@__DIR__, "lte.jl")

function parse_command_line(args)
    if isempty(args)
        return nothing
    end

    positional = String[]
    seed_mean_dir = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--seed-mean"
            i == length(args) && error("--seed-mean requires a directory path")
            seed_mean_dir = abspath(args[i + 1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end

    isempty(positional) && return nothing

    experiments_root = abspath(positional[1])
    lambda = length(positional) >= 2 ? parse(Float64, positional[2]) : 0.0
    evaluation_root = length(positional) >= 3 ? abspath(positional[3]) : joinpath(experiments_root, "super-evaluations")
    max_runs = length(positional) >= 4 ? parse(Int, positional[4]) : typemax(Int)
    return (
        experiments_root = experiments_root,
        lambda = lambda,
        evaluation_root = evaluation_root,
        max_runs = max_runs,
        seed_mean_dir = seed_mean_dir,
    )
end

function parse_ranked_sites(results_html)
    html = read(results_html, String)
    m = match(r"const ALL_SITES = (\[[\s\S]*?\]);"s, html)
    isnothing(m) && error("ALL_SITES not found in $results_html")
    sites = JSON.parse(m.captures[1])
    sort!(sites; by = s -> -Float64(s["r_val"]))
    return sites
end

function site_available(experiments_root, site)
    site_id = string(site["id"])
    site_dir = joinpath(experiments_root, site_id)
    return isdir(site_dir) &&
           isfile(joinpath(site_dir, "lt.exe.p")) &&
           isfile(joinpath(site_dir, "lte_results.csv"))
end

function available_ranked_sites(experiments_root, results_html)
    ranked = parse_ranked_sites(results_html)
    return [s for s in ranked if site_available(experiments_root, s)]
end

function scaled_cc(x)
    xf = Float64(x)
    if !isfinite(xf)
        return 1.0e-6
    end
    return clamp((xf + 1.0) / 2.0, 1.0e-6, 1.0)
end

function site_weight(site)
    return max(Float64(site["r_val"]), 0.05)
end

function weighted_mean(values, weights)
    num = sum(v * w for (v, w) in zip(values, weights))
    den = sum(weights)
    return den <= 0 ? 0.0 : num / den
end

function cumulative_score(history)
    weights = [entry["rank_weight"] for entry in history]
    val_skill = weighted_mean([scaled_cc(entry["validation_cc"]) for entry in history], weights)
    test_skill = weighted_mean([scaled_cc(entry["test_cc"]) for entry in history], weights)
    mean_coherence = weighted_mean([scaled_cc(entry["mean_manifold_cc"]) for entry in history], weights)
    penalty = weighted_mean([max(Float64(get(entry, "mean_penalty", 0.0)), 0.0) for entry in history], weights)

    # Weighted geometric mean rewards balanced performance instead of a single
    # strong metric. Validation is primary, test is secondary, manifold
    # coherence is the stability criterion across ranked indices.
    score =
        exp(0.50 * log(val_skill) +
            0.30 * log(test_skill) +
            0.20 * log(mean_coherence))

    return score - 0.05 * penalty
end

function save_json(path, value)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, value, 2)
    end
end

function copy_dir(src, dst)
    rm(dst; recursive = true, force = true)
    cp(src, dst; force = true)
end

function run_lte(experiments_root, eval_runs_root, site, lambda)
    site_id = string(site["id"])
    source_dir = joinpath(experiments_root, site_id)
    output_dir = joinpath(eval_runs_root, site_id)
    cmd = `julia --startup-file=no $(DEFAULT_LTE_SCRIPT) $source_dir $lambda $output_dir`
    run(cmd)
    summary_file = joinpath(output_dir, "optimization_summary.json")
    return JSON.parsefile(summary_file), output_dir
end

function summarize_history(history, best_entry, best_mean_dir, ranked_sites, experiments_root, eval_root; seed_mean_dir = nothing)
    return Dict(
        "experiments_root" => abspath(experiments_root),
        "evaluation_root" => abspath(eval_root),
        "seed_mean_dir" => isnothing(seed_mean_dir) ? nothing : abspath(seed_mean_dir),
        "results_html" => abspath(DEFAULT_RESULTS_HTML),
        "lte_script" => abspath(DEFAULT_LTE_SCRIPT),
        "processed_run_count" => length(history),
        "ranked_run_ids" => [string(site["id"]) for site in ranked_sites],
        "best_mean_dir" => abspath(best_mean_dir),
        "best_step" => isnothing(best_entry) ? nothing : best_entry["step"],
        "best_site_id" => isnothing(best_entry) ? nothing : best_entry["site_id"],
        "best_site_name" => isnothing(best_entry) ? nothing : best_entry["site_name"],
        "best_cumulative_score" => isnothing(best_entry) ? nothing : best_entry["cumulative_score"],
        "history" => history,
    )
end

function main()
    parsed = parse_command_line(ARGS)
    if isnothing(parsed)
        println("Usage: julia lte_super.jl <experiments_root> [lambda] [evaluation_root] [max_runs] [--seed-mean <mean_dir>]")
        return
    end

    experiments_root = parsed.experiments_root
    lambda = parsed.lambda
    evaluation_root = parsed.evaluation_root
    max_runs = parsed.max_runs
    seed_mean_dir = parsed.seed_mean_dir

    ranked_sites = available_ranked_sites(experiments_root, DEFAULT_RESULTS_HTML)
    ranked_sites = ranked_sites[1:min(length(ranked_sites), max_runs)]
    isempty(ranked_sites) && error("No ranked experiment subdirectories found under $experiments_root")

    eval_runs_root = joinpath(evaluation_root, "runs")
    snapshots_root = joinpath(evaluation_root, "snapshots")
    best_mean_dir = joinpath(evaluation_root, "best_mean")
    history_file = joinpath(evaluation_root, "super_history.json")
    summary_file = joinpath(evaluation_root, "super_summary.json")
    current_mean_dir = joinpath(eval_runs_root, "mean")

    if !isnothing(seed_mean_dir) && !isdir(seed_mean_dir)
        error("Seed mean directory not found: $seed_mean_dir")
    end
    if !isnothing(seed_mean_dir) &&
       (seed_mean_dir == evaluation_root || startswith(seed_mean_dir, evaluation_root * string(Base.Filesystem.path_separator)))
        error("Seed mean directory must be outside the new evaluation root: $seed_mean_dir")
    end

    rm(eval_runs_root; recursive = true, force = true)
    rm(snapshots_root; recursive = true, force = true)
    rm(best_mean_dir; recursive = true, force = true)
    mkpath(eval_runs_root)
    mkpath(snapshots_root)
    if !isnothing(seed_mean_dir)
        copy_dir(seed_mean_dir, current_mean_dir)
    end

    history = Vector{Dict{String, Any}}()
    best_score = -Inf
    best_entry = nothing

    println("Running ranked LTE curriculum over $(length(ranked_sites)) sites")
    println("Shared mean manifold directory: $current_mean_dir")
    if !isnothing(seed_mean_dir)
        println("Seeded initial mean manifold from: $seed_mean_dir")
    end
    flush(stdout)

    for (idx, site) in enumerate(ranked_sites)
        site_id = string(site["id"])
        site_name = string(get(site, "name", site_id))
        @printf("\n[%d/%d] %s — %s (r_val=%.3f)\n", idx, length(ranked_sites), site_id, site_name, Float64(site["r_val"]))
        flush(stdout)

        run_summary, output_dir = run_lte(experiments_root, eval_runs_root, site, lambda)
        entry = Dict{String, Any}(
            "step" => idx,
            "site_id" => site_id,
            "site_name" => site_name,
            "rank_weight" => site_weight(site),
            "seed_r_val" => Float64(site["r_val"]),
            "train_cc" => Float64(run_summary["train_cc"]),
            "validation_cc" => Float64(run_summary["validation_cc"]),
            "test_cc" => Float64(run_summary["test_cc"]),
            "mean_manifold_cc" => Float64(get(run_summary, "mean_manifold_cc", NaN)),
            "mean_penalty" => Float64(get(run_summary, "mean_penalty", 0.0)),
            "complexity" => Int(run_summary["complexity"]),
            "output_dir" => abspath(output_dir),
            "mean_dir" => abspath(current_mean_dir),
        )
        push!(history, entry)
        entry["cumulative_score"] = cumulative_score(history)

        snapshot_dir = joinpath(snapshots_root, @sprintf("%03d_%s", idx, site_id))
        copy_dir(current_mean_dir, snapshot_dir)
        entry["mean_snapshot_dir"] = abspath(snapshot_dir)

        println(
            " cumulative_score=$(round(entry["cumulative_score"], digits=6))" *
            " validation=$(round(entry["validation_cc"], digits=6))" *
            " test=$(round(entry["test_cc"], digits=6))" *
            " meanCC=$(round(entry["mean_manifold_cc"], digits=6))"
        )

        if entry["cumulative_score"] > best_score
            best_score = entry["cumulative_score"]
            best_entry = entry
            copy_dir(snapshot_dir, best_mean_dir)
            println(" -> new best mean manifold snapshot")
        end

        save_json(history_file, history)
        save_json(summary_file, summarize_history(history, best_entry, best_mean_dir, ranked_sites, experiments_root, evaluation_root; seed_mean_dir = seed_mean_dir))
    end

    println("\nBest mean manifold:")
    println(" step=$(best_entry["step"]) site=$(best_entry["site_id"]) score=$(round(best_entry["cumulative_score"], digits=6))")
    println(" saved to $best_mean_dir")
end

main()
