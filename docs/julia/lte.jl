using JSON
using DelimitedFiles
using LinearAlgebra
using Statistics

# Pukite MSL Tidal Fit Algorithm
# Ported from pukite-slr.html to Julia, with sparse evaluation support.
#
# Usage: julia lte.jl <source_dir> [lambda] [evaluation_dir]
#
# The script reads from the following files in the source directory:
# - lt.exe.p: JSON parameters (bg, delA, delB, asym, year, mp, init, lpap, ltep, harm, etc.)
# - lt.exe.resp: Response file containing YEAR adjustment
# - lte_results.csv: Data file (date, source-model, observed, manifold-target, etc.)
#
# Sparse search workflow:
# - train: outside validation and test intervals
# - validation: 1940-1970
# - test: (1970, 2000]
# - output: written to an evaluation directory so the source directory remains "gold"

const DEFAULT_IDATE = 1920.9
const VALIDATION_START = 1940.0
const VALIDATION_END = 1970.0
const TEST_START = 1970.0
const TEST_END = 2000.0
const ZERO_TOL = 1.0e-9
const DEFAULT_COORDINATE_PASSES = 2
const DEFAULT_PRUNE_PASSES = 4

mutable struct ModelState
    lpap::Matrix{Float64}
    active_lpap::BitVector
    ks::Vector{Float64}
    active_ks::BitVector
    year_len::Float64
    bg::Float64
    del_a::Float64
    del_b::Float64
    asym::Float64
    year_drift::Float64
    mp::Float64
    init_val::Float64
    ann1::Float64
    ann2::Float64
    sem1::Float64
    sem2::Float64
end

struct Dataset
    t::Vector{Float64}
    observed::Vector{Float64}
    manifold_target::Vector{Float64}
    train_mask::BitVector
    validation_mask::BitVector
    test_mask::BitVector
end

struct EvalResult
    manifold_raw::Vector{Float64}
    manifold_fit::Vector{Float64}
    annual::Vector{Float64}
    final_fit::Vector{Float64}
    beta::Vector{Float64}
    cal_intercept::Float64
    cal_slope::Float64
    train_cc::Float64
    validation_cc::Float64
    test_cc::Float64
    complexity::Int
end

"""
    build_tide(t, lpap, year_len)
Sum of 31 tidal constituents.
lpap: matrix where each row is [period_days, amplitude, phase_rad]
"""
function build_tide(t, lpap, year_len)
    N = length(t)
    tide = zeros(Float64, N)
    for c in 1:size(lpap, 1)
        pd = lpap[c, 1]
        amp = lpap[c, 2]
        phase = lpap[c, 3]
        if pd == 0 || abs(amp) <= ZERO_TOL
            continue
        end
        freq = year_len / pd
        for i in 1:N
            tide[i] += amp * cos(2π * freq * t[i] + phase)
        end
    end
    return tide
end

"""
    build_forcing(t, tide, year_drift, bg, del_a, del_b, asym)
Impulse amplification of the tidal sum.
"""
function build_forcing(t, tide, year_drift, bg, del_a, del_b, asym)
    N = length(t)
    forcing = zeros(Float64, N)
    start_t = t[1]
    d_pos = floor(Int, abs(del_b) * 12)
    h_pos = (d_pos + 6) % 12
    for i in 1:N
        dt = t[i] - start_t
        amp_time = t[i] + year_drift * dt + (bg / 1_000_000) * dt^2
        frac = amp_time - floor(amp_time)
        mi = floor(Int, frac * 12)
        imp = 0.0
        if mi == d_pos
            imp = del_a
        elseif mi == h_pos
            imp = asym
        end
        forcing[i] = tide[i] * imp
    end
    return forcing
end

"""
    build_manifold(t, forcing, idate, init_val, mp)
Bidirectional IIR integrator seeded at idate.
"""
function build_manifold(t, forcing, idate, init_val, mp)
    N = length(t)
    m = zeros(Float64, N)
    lag_c = abs(mp)

    si = findfirst(x -> x >= idate, t)
    if isnothing(si)
        si = 1
    end

    m[si] = init_val
    for i in (si + 1):N
        ramp = m[i - 1] >= 0 ? lag_c : -lag_c
        m[i] = forcing[i] + m[i - 1] - ramp
    end
    for i in si:-1:2
        ramp = m[i] >= 0 ? lag_c : -lag_c
        m[i - 1] = -forcing[i - 1] + m[i] + ramp
    end
    return m
end

"""
    build_annual(t, ann1, ann2, sem1, sem2)
Annual + semiannual fitted cycles.
"""
function build_annual(t, ann1, ann2, sem1, sem2)
    N = length(t)
    ann = zeros(Float64, N)
    for i in 1:N
        ann[i] = ann1 * cos(2π * t[i] + ann2) + sem1 * cos(4π * t[i] + sem2)
    end
    return ann
end

function build_design_row(mi, ks)
    row = [1.0, mi]
    for k in ks
        push!(row, sin(2π * k * mi))
        push!(row, cos(2π * k * mi))
    end
    return row
end

"""
    solve_mlr(manifold, y_resid, ks, lambda)
Multiple Linear Regression using normal equations with Ridge regularization.
"""
function solve_mlr(manifold, y_resid, ks, lambda)
    N = length(manifold)
    P = 2 + length(ks) * 2
    XtX = zeros(Float64, P, P)
    Xty = zeros(Float64, P)

    for i in 1:N
        row = build_design_row(manifold[i], ks)
        yi = y_resid[i]
        for j in 1:P
            Xty[j] += row[j] * yi
            for k in 1:P
                XtX[j, k] += row[j] * row[k]
            end
        end
    end

    for j in 2:P
        XtX[j, j] += lambda
    end

    return XtX \ Xty
end

function apply_beta(manifold, beta, ks)
    N = length(manifold)
    fit = zeros(Float64, N)
    for i in 1:N
        row = build_design_row(manifold[i], ks)
        s = 0.0
        for j in 1:length(beta)
            s += beta[j] * row[j]
        end
        fit[i] = s
    end
    return fit
end

function solve_linear_map(x_vals, y_vals)
    n = length(x_vals)
    sx = sum(x_vals)
    sy = sum(y_vals)
    sxx = sum(x_vals .^ 2)
    sxy = sum(x_vals .* y_vals)

    denom = n * sxx - sx^2
    if abs(denom) < 1e-15
        return (intercept = sy / n, slope = 0.0)
    end
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return (intercept = intercept, slope = slope)
end

function apply_linear_map(x_vals, map)
    return map.intercept .+ map.slope .* x_vals
end

function safe_corrcoef(x_vals, y_vals)
    if length(x_vals) < 2 || length(y_vals) < 2
        return NaN
    end
    sx = std(x_vals)
    sy = std(y_vals)
    if !isfinite(sx) || !isfinite(sy) || sx <= ZERO_TOL || sy <= ZERO_TOL
        return NaN
    end
    return cor(x_vals, y_vals)
end

function parse_resp_year(resp_file)
    year_adj = 0.005
    if isfile(resp_file)
        for line in readlines(resp_file)
            parts = split(line)
            if length(parts) >= 2 && parts[1] == "YEAR"
                val = tryparse(Float64, replace(parts[2], "\"" => ""))
                if !isnothing(val)
                    year_adj = val
                end
            end
        end
    end
    return year_adj
end

function derive_ks(params)
    if haskey(params, "sparse") && params["sparse"] isa AbstractDict && haskey(params["sparse"], "ks")
        return Float64.(params["sparse"]["ks"])
    end

    ltep = Float64.(params["ltep"])
    harm = Float64.(params["harm"])
    ks = Float64[]
    if !isempty(ltep)
        push!(ks, ltep[1])
    end
    if length(ltep) >= 2
        push!(ks, ltep[2])
        for h in harm
            push!(ks, ltep[2] * Float64(h))
        end
    end
    return ks
end

function load_lpap_matrix(params)
    lpap_raw = params["lpap"]
    lpap_mat = zeros(Float64, length(lpap_raw), 3)
    for i in 1:length(lpap_raw)
        lpap_mat[i, :] = [Float64(x) for x in lpap_raw[i]]
    end
    return lpap_mat
end

function initial_state(params, year_adj)
    lpap = load_lpap_matrix(params)
    active_lpap = BitVector(abs.(lpap[:, 2]) .> ZERO_TOL)
    if haskey(params, "sparse") && params["sparse"] isa AbstractDict && haskey(params["sparse"], "active_lpap")
        active_lpap = BitVector(Bool.(params["sparse"]["active_lpap"]))
    end

    ks = derive_ks(params)
    active_ks = trues(length(ks))
    if haskey(params, "sparse") && params["sparse"] isa AbstractDict && haskey(params["sparse"], "active_ks")
        active_ks = BitVector(Bool.(params["sparse"]["active_ks"]))
    end

    return ModelState(
        lpap,
        active_lpap,
        ks,
        active_ks,
        365.2422484 + year_adj,
        Float64(params["bg"]),
        Float64(params["delA"]),
        Float64(params["delB"]),
        Float64(params["asym"]),
        Float64(params["year"]),
        Float64(params["mp"]),
        Float64(params["init"]),
        Float64(params["ann1"]),
        Float64(params["ann2"]),
        Float64(params["sem1"]),
        Float64(params["sem2"]),
    )
end

function copy_state(state)
    return ModelState(
        copy(state.lpap),
        copy(state.active_lpap),
        copy(state.ks),
        copy(state.active_ks),
        state.year_len,
        state.bg,
        state.del_a,
        state.del_b,
        state.asym,
        state.year_drift,
        state.mp,
        state.init_val,
        state.ann1,
        state.ann2,
        state.sem1,
        state.sem2,
    )
end

function active_lpap_matrix(state)
    lpap = copy(state.lpap)
    for i in 1:size(lpap, 1)
        if !state.active_lpap[i]
            lpap[i, 2] = 0.0
        end
    end
    return lpap
end

function active_ks(state)
    return [state.ks[i] for i in eachindex(state.ks) if state.active_ks[i]]
end

function make_masks(t)
    validation_mask = BitVector((t .>= VALIDATION_START) .& (t .<= VALIDATION_END))
    test_mask = BitVector((t .> TEST_START) .& (t .<= TEST_END))
    train_mask = .!(validation_mask .| test_mask)
    return train_mask, validation_mask, test_mask
end

function load_dataset(csv_file)
    data = readdlm(csv_file, ',')
    t = Float64.(data[:, 1])
    observed = Float64.(data[:, 3])
    manifold_target = Float64.(data[:, 4])
    train_mask, validation_mask, test_mask = make_masks(t)
    return Dataset(t, observed, manifold_target, train_mask, validation_mask, test_mask)
end

function count_complexity(state)
    complexity = count(identity, state.active_lpap) + count(identity, state.active_ks)
    complexity += abs(state.bg) > ZERO_TOL ? 1 : 0
    complexity += abs(state.del_a) > ZERO_TOL ? 1 : 0
    complexity += abs(state.asym) > ZERO_TOL ? 1 : 0
    complexity += abs(state.year_drift) > ZERO_TOL ? 1 : 0
    complexity += abs(state.mp) > ZERO_TOL ? 1 : 0
    complexity += abs(state.init_val) > ZERO_TOL ? 1 : 0
    complexity += abs(state.ann1) > ZERO_TOL ? 1 : 0
    complexity += abs(state.sem1) > ZERO_TOL ? 1 : 0
    return complexity
end

function evaluate_state(dataset, state, lambda)
    lpap = active_lpap_matrix(state)
    ks = active_ks(state)
    tide = build_tide(dataset.t, lpap, state.year_len)
    forcing = build_forcing(dataset.t, tide, state.year_drift, state.bg, state.del_a, state.del_b, state.asym)
    manifold_raw = build_manifold(dataset.t, forcing, DEFAULT_IDATE, state.init_val, state.mp)
    annual = build_annual(dataset.t, state.ann1, state.ann2, state.sem1, state.sem2)

    cal = solve_linear_map(manifold_raw, dataset.manifold_target)
    manifold_fit = apply_linear_map(manifold_raw, cal)

    y_resid = dataset.observed .- annual
    beta = solve_mlr(manifold_fit[dataset.train_mask], y_resid[dataset.train_mask], ks, lambda)
    fit_resid = apply_beta(manifold_fit, beta, ks)
    final_fit = fit_resid .+ annual

    return EvalResult(
        manifold_raw,
        manifold_fit,
        annual,
        final_fit,
        beta,
        cal.intercept,
        cal.slope,
        safe_corrcoef(dataset.observed[dataset.train_mask], final_fit[dataset.train_mask]),
        safe_corrcoef(dataset.observed[dataset.validation_mask], final_fit[dataset.validation_mask]),
        safe_corrcoef(dataset.observed[dataset.test_mask], final_fit[dataset.test_mask]),
        count_complexity(state),
    )
end

score_value(x) = isfinite(x) ? x : -1.0

function better_result(candidate, incumbent; tol = 1.0e-6)
    if score_value(candidate.validation_cc) > score_value(incumbent.validation_cc) + tol
        return true
    elseif score_value(candidate.validation_cc) + tol < score_value(incumbent.validation_cc)
        return false
    end

    if score_value(candidate.train_cc) > score_value(incumbent.train_cc) + tol
        return true
    elseif score_value(candidate.train_cc) + tol < score_value(incumbent.train_cc)
        return false
    end

    return candidate.complexity < incumbent.complexity
end

function param_step(name, value)
    if name == :del_b
        return max(abs(value) * 0.05, 0.01)
    elseif name in (:ann2, :sem2)
        return 0.1
    elseif name in (:bg, :mp, :init_val, :year_drift)
        return max(abs(value) * 0.25, 1.0e-5)
    else
        return max(abs(value) * 0.15, 1.0e-3)
    end
end

function coordinate_specs(state)
    specs = Any[]
    for i in eachindex(state.active_lpap)
        if state.active_lpap[i]
            push!(specs, (kind = :lpap_amp, index = i))
            push!(specs, (kind = :lpap_phase, index = i))
        end
    end
    for i in eachindex(state.ks)
        if state.active_ks[i]
            push!(specs, (kind = :ks, index = i))
        end
    end
    for name in (:bg, :del_a, :del_b, :asym, :year_drift, :mp, :init_val, :ann1, :ann2, :sem1, :sem2)
        push!(specs, (kind = :param, name = name))
    end
    return specs
end

function current_value(state, spec)
    if spec.kind == :lpap_amp
        return state.lpap[spec.index, 2]
    elseif spec.kind == :lpap_phase
        return state.lpap[spec.index, 3]
    elseif spec.kind == :ks
        return state.ks[spec.index]
    else
        return getfield(state, spec.name)
    end
end

function current_step(state, spec)
    if spec.kind == :lpap_amp
        return max(abs(state.lpap[spec.index, 2]) * 0.15, 1.0e-4)
    elseif spec.kind == :lpap_phase
        return 0.1
    elseif spec.kind == :ks
        return max(abs(state.ks[spec.index]) * 0.05, 0.05)
    else
        return param_step(spec.name, getfield(state, spec.name))
    end
end

function set_value!(state, spec, value)
    if spec.kind == :lpap_amp
        state.lpap[spec.index, 2] = value
        state.active_lpap[spec.index] = abs(value) > ZERO_TOL
    elseif spec.kind == :lpap_phase
        state.lpap[spec.index, 3] = value
    elseif spec.kind == :ks
        state.ks[spec.index] = value
        state.active_ks[spec.index] = abs(value) > ZERO_TOL
    else
        setfield!(state, spec.name, value)
    end
end

function coordinate_refine(dataset, state, result, lambda; passes = DEFAULT_COORDINATE_PASSES)
    best_state = copy_state(state)
    best_result = result
    for _ in 1:passes
        improved = false
        for spec in coordinate_specs(best_state)
            base_value = current_value(best_state, spec)
            step = current_step(best_state, spec)
            candidates = (base_value + step, base_value - step)
            for new_value in candidates
                if spec.kind == :ks && new_value <= ZERO_TOL
                    continue
                end
                candidate_state = copy_state(best_state)
                set_value!(candidate_state, spec, new_value)
                candidate_result = evaluate_state(dataset, candidate_state, lambda)
                if better_result(candidate_result, best_result)
                    best_state = candidate_state
                    best_result = candidate_result
                    improved = true
                    break
                end
            end
        end
        if !improved
            break
        end
    end
    return best_state, best_result
end

function prune_targets(state)
    targets = Any[]
    for i in eachindex(state.active_lpap)
        if state.active_lpap[i]
            push!(targets, (kind = :lpap, index = i))
        end
    end
    for i in eachindex(state.active_ks)
        if state.active_ks[i]
            push!(targets, (kind = :ks, index = i))
        end
    end
    for name in (:bg, :del_a, :asym, :year_drift, :mp, :init_val, :ann1, :sem1)
        if abs(getfield(state, name)) > ZERO_TOL
            push!(targets, (kind = :param, name = name))
        end
    end
    return targets
end

function apply_prune!(state, target)
    if target.kind == :lpap
        state.active_lpap[target.index] = false
        state.lpap[target.index, 2] = 0.0
    elseif target.kind == :ks
        state.active_ks[target.index] = false
    else
        setfield!(state, target.name, 0.0)
        if target.name == :ann1
            state.ann2 = 0.0
        elseif target.name == :sem1
            state.sem2 = 0.0
        end
    end
end

function sparse_optimize(dataset, state, lambda)
    best_state = copy_state(state)
    best_result = evaluate_state(dataset, best_state, lambda)
    best_state, best_result = coordinate_refine(dataset, best_state, best_result, lambda)

    for _ in 1:DEFAULT_PRUNE_PASSES
        trial_state = nothing
        trial_result = nothing
        for target in prune_targets(best_state)
            candidate_state = copy_state(best_state)
            apply_prune!(candidate_state, target)
            candidate_result = evaluate_state(dataset, candidate_state, lambda)
            candidate_state, candidate_result = coordinate_refine(dataset, candidate_state, candidate_result, lambda; passes = 1)
            if isnothing(trial_result) || better_result(candidate_result, trial_result)
                trial_state = candidate_state
                trial_result = candidate_result
            end
        end

        if isnothing(trial_result) || !better_result(trial_result, best_result)
            break
        end

        best_state = trial_state
        best_result = trial_result
    end

    return best_state, best_result
end

function default_evaluation_dir(source_dir)
    abs_source = abspath(source_dir)
    if basename(dirname(abs_source)) == "evaluations"
        return abs_source
    end
    return joinpath(dirname(abs_source), "evaluations", basename(abs_source))
end

function state_to_params(params, state, source_dir, result)
    out = deepcopy(params)
    masked_lpap = active_lpap_matrix(state)
    out["lpap"] = [[masked_lpap[i, 1], masked_lpap[i, 2], masked_lpap[i, 3]] for i in 1:size(masked_lpap, 1)]
    out["bg"] = state.bg
    out["delA"] = state.del_a
    out["delB"] = state.del_b
    out["asym"] = state.asym
    out["year"] = state.year_drift
    out["mp"] = state.mp
    out["init"] = state.init_val
    out["ann1"] = state.ann1
    out["ann2"] = state.ann2
    out["sem1"] = state.sem1
    out["sem2"] = state.sem2
    out["sparse"] = Dict(
        "source_dir" => abspath(source_dir),
        "validation_interval" => [VALIDATION_START, VALIDATION_END],
        "test_interval" => [TEST_START, TEST_END],
        "ks" => active_ks(state),
        "active_lpap" => collect(state.active_lpap),
        "active_ks" => collect(trues(length(active_ks(state)))),
        "method" => "greedy-pruning-plus-coordinate-refinement",
        "objective" => "validation_cc primary, training_cc tie-breaker",
        "metrics" => Dict(
            "train_cc" => result.train_cc,
            "validation_cc" => result.validation_cc,
            "test_cc" => result.test_cc,
            "complexity" => result.complexity,
        ),
    )
    return out
end

function save_outputs(source_dir, output_dir, params, state, result, dataset)
    mkpath(output_dir)

    src_csv = joinpath(source_dir, "lte_results.csv")
    src_resp = joinpath(source_dir, "lt.exe.resp")
    dst_csv = joinpath(output_dir, "lte_results.csv")
    dst_resp = joinpath(output_dir, "lt.exe.resp")

    if abspath(src_csv) != abspath(dst_csv)
        cp(src_csv, dst_csv; force = true)
    end
    if isfile(src_resp) && abspath(src_resp) != abspath(dst_resp)
        cp(src_resp, dst_resp; force = true)
    end

    out_params = state_to_params(params, state, source_dir, result)
    open(joinpath(output_dir, "lt.exe.p"), "w") do io
        JSON.print(io, out_params, 2)
    end

    results_matrix = hcat(dataset.t, result.manifold_raw, result.manifold_fit, result.final_fit)
    writedlm(joinpath(output_dir, "julia_results.csv"), results_matrix, ',')

    summary = Dict(
        "source_dir" => abspath(source_dir),
        "evaluation_dir" => abspath(output_dir),
        "validation_interval" => [VALIDATION_START, VALIDATION_END],
        "test_interval" => [TEST_START, TEST_END],
        "train_cc" => result.train_cc,
        "validation_cc" => result.validation_cc,
        "test_cc" => result.test_cc,
        "complexity" => result.complexity,
        "active_lpap_count" => count(identity, state.active_lpap),
        "active_ks" => active_ks(state),
        "beta" => result.beta,
        "manifold_calibration" => Dict(
            "intercept" => result.cal_intercept,
            "slope" => result.cal_slope,
        ),
    )
    open(joinpath(output_dir, "optimization_summary.json"), "w") do io
        JSON.print(io, summary, 2)
    end
end

function main()
    if length(ARGS) < 1
        println("Usage: julia lte.jl <source_dir> [lambda] [evaluation_dir]")
        return
    end

    source_dir = abspath(ARGS[1])
    lambda = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0
    output_dir = length(ARGS) >= 3 ? abspath(ARGS[3]) : default_evaluation_dir(source_dir)

    p_file = joinpath(source_dir, "lt.exe.p")
    resp_file = joinpath(source_dir, "lt.exe.resp")
    csv_file = joinpath(source_dir, "lte_results.csv")

    if !isfile(p_file) || !isfile(csv_file)
        println("Error: Required files not found in $source_dir")
        return
    end

    params = JSON.parsefile(p_file)
    year_adj = parse_resp_year(resp_file)
    dataset = load_dataset(csv_file)
    state = initial_state(params, year_adj)
    best_state, best_result = sparse_optimize(dataset, state, lambda)
    save_outputs(source_dir, output_dir, params, best_state, best_result, dataset)

    println(
        "CC train=$(round(best_result.train_cc, digits=6)) " *
        "validation=$(round(best_result.validation_cc, digits=6)) " *
        "test=$(round(best_result.test_cc, digits=6)) " *
        "complexity=$(best_result.complexity)"
    )
    println("Sparse evaluation saved to $output_dir")
end

main()
