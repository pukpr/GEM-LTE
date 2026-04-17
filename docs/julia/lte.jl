using JSON
using DelimitedFiles
using LinearAlgebra
using Statistics

# Pukite MSL Tidal Fit Algorithm
# Ported from pukite-slr.html to Julia, with sparse evaluation support.
#
# Usage: julia lte.jl <source_dir> [lambda] [evaluation_dir] [--seed-mean <mean_dir>] [--strict]
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
# - shared mean manifold: ../evaluations/mean accumulates the averaged LPAP amplitudes across runs

const DEFAULT_IDATE = 1920.9
const VALIDATION_START = 1970.0
const VALIDATION_END = 2000.0
const TEST_START = 1940.0
const TEST_END = 1970.0
const ZERO_TOL = 1.0e-9
const GAP_TOL = 1.0e-4
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
    observed_mask::BitVector
    target_mask::BitVector
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
    mean_manifold_cc::Float64
    strict_mean_manifold_cc::Float64
    mean_penalty_weight::Float64
    mean_deviation::Float64
    mean_penalty::Float64
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
    X = zeros(Float64, N, P)

    for i in 1:N
        row = build_design_row(manifold[i], ks)
        for j in 1:P
            X[i, j] = row[j]
        end
    end

    if lambda > 0.0
        X_aug = zeros(Float64, N + P, P)
        y_aug = zeros(Float64, N + P)
        X_aug[1:N, :] .= X
        y_aug[1:N] .= y_resid
        for j in 2:P
            X_aug[N + j, j] = sqrt(lambda)
        end
        return X_aug \ y_aug
    end

    return X \ y_resid
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

function safe_corrcoef_masked(x_vals, y_vals, mask)
    idx = findall(mask .& (abs.(x_vals) .> ZERO_TOL) .& (abs.(y_vals) .> ZERO_TOL))
    if length(idx) < 2
        return NaN
    end
    return safe_corrcoef(x_vals[idx], y_vals[idx])
end

function infer_month_step(t)
    diffs = [t[i + 1] - t[i] for i in 1:(length(t) - 1) if t[i + 1] > t[i]]
    if isempty(diffs)
        return 1.0 / 12.0
    end
    return median(diffs)
end

function expand_monthly_grid(t, observed, manifold_target)
    if length(t) <= 1
        observed_mask = BitVector(observed .!= 0.0)
        target_mask = trues(length(t))
        return t, observed, manifold_target, observed_mask, target_mask
    end

    step = infer_month_step(t)
    expanded_t = Float64[t[1]]
    expanded_observed = Float64[observed[1]]
    expanded_target = Float64[manifold_target[1]]
    observed_mask = BitVector([observed[1] != 0.0])
    target_mask = BitVector([true])

    for i in 2:length(t)
        prev_t = t[i - 1]
        dt = t[i] - prev_t
        if dt > step + GAP_TOL
            missing_n = max(round(Int, dt / step) - 1, 0)
            for j in 1:missing_n
                push!(expanded_t, prev_t + step * j)
                push!(expanded_observed, 0.0)
                push!(expanded_target, 0.0)
                push!(observed_mask, false)
                push!(target_mask, false)
            end
        end
        push!(expanded_t, t[i])
        push!(expanded_observed, observed[i])
        push!(expanded_target, manifold_target[i])
        push!(observed_mask, observed[i] != 0.0)
        push!(target_mask, true)
    end

    return expanded_t, expanded_observed, expanded_target, observed_mask, target_mask
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

function resolve_gold_source_dir(source_dir, params)
    if haskey(params, "sparse") && params["sparse"] isa AbstractDict && haskey(params["sparse"], "source_dir")
        sparse_source = String(params["sparse"]["source_dir"])
        if isdir(sparse_source)
            return abspath(sparse_source)
        end
    end
    return abspath(source_dir)
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

function load_state_from_dir(source_dir)
    p_file = joinpath(source_dir, "lt.exe.p")
    resp_file = joinpath(source_dir, "lt.exe.resp")
    params = JSON.parsefile(p_file)
    year_adj = parse_resp_year(resp_file)
    return params, initial_state(params, year_adj)
end

function state_with_mean_manifold(run_state, mean_state)
    merged = copy_state(run_state)
    n_lpap = min(size(merged.lpap, 1), size(mean_state.lpap, 1))
    for i in 1:n_lpap
        merged.lpap[i, 2] = abs(mean_state.lpap[i, 2]) > ZERO_TOL ? mean_state.lpap[i, 2] : 0.0
        merged.active_lpap[i] = abs(merged.lpap[i, 2]) > ZERO_TOL
    end
    return merged
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

function lpap_amplitude_vector(state)
    amps = zeros(Float64, size(state.lpap, 1))
    for i in 1:length(amps)
        if state.active_lpap[i]
            amps[i] = state.lpap[i, 2]
        end
    end
    return amps
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
    t_raw = Float64.(data[:, 1])
    observed_raw = Float64.(data[:, 3])
    manifold_target_raw = Float64.(data[:, 4])
    t, observed, manifold_target, observed_mask, target_mask =
        expand_monthly_grid(t_raw, observed_raw, manifold_target_raw)
    train_mask, validation_mask, test_mask = make_masks(t)
    return Dataset(t, observed, manifold_target, observed_mask, target_mask, train_mask, validation_mask, test_mask)
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

function manifold_raw_for_state(t, state)
    lpap = active_lpap_matrix(state)
    tide = build_tide(t, lpap, state.year_len)
    forcing = build_forcing(t, tide, state.year_drift, state.bg, state.del_a, state.del_b, state.asym)
    return build_manifold(t, forcing, DEFAULT_IDATE, state.init_val, state.mp)
end

function strict_mean_corr(manifold_raw, t, mean_state)
    return safe_corrcoef(manifold_raw, manifold_raw_for_state(t, mean_state))
end

function bar_chart_corr(state, mean_state)
    return safe_corrcoef(lpap_amplitude_vector(state), lpap_amplitude_vector(mean_state))
end

function training_mask(dataset)
    return dataset.train_mask .& dataset.observed_mask
end

function validation_mask(dataset)
    return dataset.validation_mask .& dataset.observed_mask
end

function test_mask(dataset)
    return dataset.test_mask .& dataset.observed_mask
end

phase_distance(a, b) = abs(atan(sin(a - b), cos(a - b)))

function mean_penalty_weight(mean_run_count)
    return mean_run_count <= 1 ? 0.0 : min(0.20, 0.02 * (mean_run_count - 1))
end

function strict_manifold_deviation_from_mean(state, mean_state)
    terms = Float64[]
    current_lpap = active_lpap_matrix(state)
    mean_lpap = active_lpap_matrix(mean_state)
    n_lpap = min(size(current_lpap, 1), size(mean_lpap, 1))
    for i in 1:n_lpap
        amp_scale = max(abs(mean_lpap[i, 2]), 1.0e-3)
        push!(terms, abs(current_lpap[i, 2] - mean_lpap[i, 2]) / amp_scale)
        push!(terms, phase_distance(current_lpap[i, 3], mean_lpap[i, 3]) / π)
    end
    for field in (:year_len, :bg, :del_a, :del_b, :asym, :year_drift, :mp, :init_val, :ann1, :ann2, :sem1, :sem2)
        mean_val = getfield(mean_state, field)
        cur_val = getfield(state, field)
        scale = field in (:ann2, :sem2) ? π : max(abs(mean_val), 1.0e-3)
        push!(terms, abs(cur_val - mean_val) / scale)
    end
    return isempty(terms) ? 0.0 : mean(terms)
end

function bar_chart_deviation_from_mean(state, mean_state)
    terms = Float64[]
    current_amps = lpap_amplitude_vector(state)
    mean_amps = lpap_amplitude_vector(mean_state)
    n_lpap = min(length(current_amps), length(mean_amps))
    for i in 1:n_lpap
        amp_scale = max(abs(mean_amps[i]), 1.0e-3)
        push!(terms, abs(current_amps[i] - mean_amps[i]) / amp_scale)
    end
    return isempty(terms) ? 0.0 : mean(terms)
end

function effective_validation_score(result)
    return score_value(result.validation_cc) - result.mean_penalty
end

function effective_training_score(result)
    return score_value(result.train_cc) - 0.5 * result.mean_penalty
end

function evaluate_state(dataset, state, lambda; mean_state = nothing, mean_run_count = 0, strict = false)
    ks = active_ks(state)
    manifold_raw = manifold_raw_for_state(dataset.t, state)
    annual = build_annual(dataset.t, state.ann1, state.ann2, state.sem1, state.sem2)

    cal_idx = findall(dataset.target_mask)
    cal = solve_linear_map(manifold_raw[cal_idx], dataset.manifold_target[cal_idx])
    manifold_fit = apply_linear_map(manifold_raw, cal)

    y_resid = dataset.observed .- annual
    train_idx = findall(training_mask(dataset))
    beta = solve_mlr(manifold_fit[train_idx], y_resid[train_idx], ks, lambda)
    fit_resid = apply_beta(manifold_fit, beta, ks)
    final_fit = fit_resid .+ annual

    mean_manifold_cc = NaN
    strict_mean_manifold_cc = NaN
    deviation = 0.0
    penalty_weight = 0.0
    if !isnothing(mean_state)
        strict_mean_manifold_cc = strict_mean_corr(manifold_raw, dataset.t, mean_state)
        mean_manifold_cc = strict ? strict_mean_manifold_cc : bar_chart_corr(state, mean_state)
        deviation = strict ? strict_manifold_deviation_from_mean(state, mean_state) : bar_chart_deviation_from_mean(state, mean_state)
        penalty_weight = mean_penalty_weight(mean_run_count)
    end
    penalty = penalty_weight * deviation

    return EvalResult(
        manifold_raw,
        manifold_fit,
        annual,
        final_fit,
        beta,
        cal.intercept,
        cal.slope,
        safe_corrcoef_masked(dataset.observed, final_fit, training_mask(dataset)),
        safe_corrcoef_masked(dataset.observed, final_fit, validation_mask(dataset)),
        safe_corrcoef_masked(dataset.observed, final_fit, test_mask(dataset)),
        mean_manifold_cc,
        strict_mean_manifold_cc,
        penalty_weight,
        deviation,
        penalty,
        count_complexity(state),
    )
end

score_value(x) = isfinite(x) ? x : -1.0

function better_result(candidate, incumbent; tol = 1.0e-6)
    if effective_validation_score(candidate) > effective_validation_score(incumbent) + tol
        return true
    elseif effective_validation_score(candidate) + tol < effective_validation_score(incumbent)
        return false
    end

    if effective_training_score(candidate) > effective_training_score(incumbent) + tol
        return true
    elseif effective_training_score(candidate) + tol < effective_training_score(incumbent)
        return false
    end

    if score_value(candidate.mean_manifold_cc) > score_value(incumbent.mean_manifold_cc) + tol
        return true
    elseif score_value(candidate.mean_manifold_cc) + tol < score_value(incumbent.mean_manifold_cc)
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

function coordinate_refine(dataset, state, result, lambda; passes = DEFAULT_COORDINATE_PASSES, mean_state = nothing, mean_run_count = 0, strict = false)
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
                candidate_result = evaluate_state(dataset, candidate_state, lambda; mean_state = mean_state, mean_run_count = mean_run_count, strict = strict)
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

function sparse_optimize(dataset, state, lambda; mean_state = nothing, mean_run_count = 0, strict = false)
    best_state = copy_state(state)
    best_result = evaluate_state(dataset, best_state, lambda; mean_state = mean_state, mean_run_count = mean_run_count, strict = strict)
    best_state, best_result = coordinate_refine(dataset, best_state, best_result, lambda; mean_state = mean_state, mean_run_count = mean_run_count, strict = strict)

    for _ in 1:DEFAULT_PRUNE_PASSES
        trial_state = nothing
        trial_result = nothing
        for target in prune_targets(best_state)
            candidate_state = copy_state(best_state)
            apply_prune!(candidate_state, target)
            candidate_result = evaluate_state(dataset, candidate_state, lambda; mean_state = mean_state, mean_run_count = mean_run_count, strict = strict)
            candidate_state, candidate_result = coordinate_refine(dataset, candidate_state, candidate_result, lambda; passes = 1, mean_state = mean_state, mean_run_count = mean_run_count, strict = strict)
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

function mean_evaluation_dir(output_dir)
    return joinpath(dirname(abspath(output_dir)), "mean")
end

function mean_run_count(mean_dir)
    summary_file = joinpath(mean_dir, "optimization_summary.json")
    if isfile(summary_file)
        summary = JSON.parsefile(summary_file)
        if haskey(summary, "mean_run_count")
            return Int(summary["mean_run_count"])
        end
    end
    return 0
end

function mean_state_from_dir(mean_dir)
    p_file = joinpath(mean_dir, "lt.exe.p")
    if !isfile(p_file)
        return nothing, 0
    end
    _, state = load_state_from_dir(mean_dir)
    return state, mean_run_count(mean_dir)
end

function averaged_mean_state(old_mean_state, old_count, current_state)
    if old_count <= 0 || isnothing(old_mean_state)
        return copy_state(current_state), 1
    end

    total = old_count + 1
    mean_state = copy_state(old_mean_state)
    current_amps = lpap_amplitude_vector(current_state)
    old_amps = lpap_amplitude_vector(old_mean_state)

    for i in 1:size(mean_state.lpap, 1)
        mean_state.lpap[i, 1] = old_mean_state.lpap[i, 1]
        mean_state.lpap[i, 2] = (old_amps[i] * old_count + current_amps[i]) / total
        mean_state.active_lpap[i] = abs(mean_state.lpap[i, 2]) > ZERO_TOL
    end
    return mean_state, total
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

function save_outputs(source_dir, output_dir, params, state, result, dataset; mean_dir = nothing, mean_run_count = 0, strict = false)
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
        "mean_dir" => isnothing(mean_dir) ? nothing : abspath(mean_dir),
        "mean_run_count" => mean_run_count,
        "mean_comparison_mode" => strict ? "strict-manifold" : "lpap-bar-chart",
        "validation_interval" => [VALIDATION_START, VALIDATION_END],
        "test_interval" => [TEST_START, TEST_END],
        "train_cc" => result.train_cc,
        "validation_cc" => result.validation_cc,
        "test_cc" => result.test_cc,
        "mean_manifold_cc" => result.mean_manifold_cc,
        "strict_mean_manifold_cc" => result.strict_mean_manifold_cc,
        "mean_penalty_weight" => result.mean_penalty_weight,
        "mean_deviation" => result.mean_deviation,
        "mean_penalty" => result.mean_penalty,
        "complexity" => result.complexity,
        "observed_count" => count(identity, dataset.observed_mask),
        "missing_observation_count" => length(dataset.observed_mask) - count(identity, dataset.observed_mask),
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

function save_mean_outputs(mean_dir, template_params, source_dir, state, run_count)
    mkpath(mean_dir)
    params = deepcopy(template_params)
    masked_lpap = active_lpap_matrix(state)
    params["lpap"] = [[masked_lpap[i, 1], masked_lpap[i, 2], masked_lpap[i, 3]] for i in 1:size(masked_lpap, 1)]
    params["bg"] = state.bg
    params["delA"] = state.del_a
    params["delB"] = state.del_b
    params["asym"] = state.asym
    params["year"] = state.year_drift
    params["mp"] = state.mp
    params["init"] = state.init_val
    params["ann1"] = state.ann1
    params["ann2"] = state.ann2
    params["sem1"] = state.sem1
    params["sem2"] = state.sem2
    params["sparse"] = Dict(
        "source_dir" => abspath(source_dir),
        "mean_run_count" => run_count,
        "mean_role" => "shared-lpap-amplitudes",
        "ks" => active_ks(state),
        "active_lpap" => collect(state.active_lpap),
        "active_ks" => collect(state.active_ks),
    )

    open(joinpath(mean_dir, "lt.exe.p"), "w") do io
        JSON.print(io, params, 2)
    end
    open(joinpath(mean_dir, "lt.exe.resp"), "w") do io
        year_adj = state.year_len - 365.2422484
        write(io, "YEAR $(year_adj)\n")
    end
    summary = Dict(
        "mean_dir" => abspath(mean_dir),
        "source_dir" => abspath(source_dir),
        "mean_run_count" => run_count,
        "active_lpap_count" => count(identity, state.active_lpap),
        "active_ks" => active_ks(state),
    )
    open(joinpath(mean_dir, "optimization_summary.json"), "w") do io
        JSON.print(io, summary, 2)
    end
end

function fmt_num(x; digits = 6)
    return string(round(x, digits = digits))
end

function pruned_summary(gold_state, best_state)
    pruned_lpap = String[]
    for i in eachindex(gold_state.active_lpap)
        if gold_state.active_lpap[i] && (i > length(best_state.active_lpap) || !best_state.active_lpap[i])
            push!(pruned_lpap, "period=$(fmt_num(gold_state.lpap[i, 1]))d")
        end
    end

    pruned_ks = String[]
    for i in eachindex(gold_state.active_ks)
        if gold_state.active_ks[i] && (i > length(best_state.active_ks) || !best_state.active_ks[i])
            push!(pruned_ks, fmt_num(gold_state.ks[i]))
        end
    end

    param_specs = [
        (:bg, "bg"),
        (:del_a, "delA"),
        (:asym, "asym"),
        (:year_drift, "year"),
        (:mp, "mp"),
        (:init_val, "init"),
        (:ann1, "ann1"),
        (:sem1, "sem1"),
    ]
    pruned_params = String[]
    for (field, label) in param_specs
        if abs(getfield(gold_state, field)) > ZERO_TOL && abs(getfield(best_state, field)) <= ZERO_TOL
            push!(pruned_params, label)
        end
    end

    lpap_text = isempty(pruned_lpap) ? "none" : join(pruned_lpap, ", ")
    ks_text = isempty(pruned_ks) ? "none" : join(pruned_ks, ", ")
    params_text = isempty(pruned_params) ? "none" : join(pruned_params, ", ")
    return lpap_text, ks_text, params_text
end

function beta_summary(result, state)
    ks = active_ks(state)
    parts = ["β0(intercept)=$(fmt_num(result.beta[1]))", "βm(manifold)=$(fmt_num(result.beta[2]))"]
    for (j, k) in enumerate(ks)
        sin_idx = 2 + 2 * j - 1
        cos_idx = sin_idx + 1
        push!(
            parts,
            "k=$(fmt_num(k)) => βsin=$(fmt_num(result.beta[sin_idx])) βcos=$(fmt_num(result.beta[cos_idx]))"
        )
    end
    return join(parts, "; ")
end

function parse_command_line(args)
    if isempty(args)
        return nothing
    end

    positional = String[]
    seed_mean_dir = nothing
    strict = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--seed-mean"
            i == length(args) && error("--seed-mean requires a directory path")
            seed_mean_dir = abspath(args[i + 1])
            i += 2
        elseif arg == "--strict"
            strict = true
            i += 1
        else
            push!(positional, arg)
            i += 1
        end
    end

    isempty(positional) && return nothing

    source_dir = abspath(positional[1])
    lambda = length(positional) >= 2 ? parse(Float64, positional[2]) : 0.0
    output_dir = length(positional) >= 3 ? abspath(positional[3]) : default_evaluation_dir(source_dir)
    return (source_dir = source_dir, lambda = lambda, output_dir = output_dir, seed_mean_dir = seed_mean_dir, strict = strict)
end

function main()
    parsed = parse_command_line(ARGS)
    if isnothing(parsed)
        println("Usage: julia lte.jl <source_dir> [lambda] [evaluation_dir] [--seed-mean <mean_dir>] [--strict]")
        return
    end

    source_dir = parsed.source_dir
    lambda = parsed.lambda
    output_dir = parsed.output_dir
    seed_mean_dir = parsed.seed_mean_dir
    strict = parsed.strict

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
    gold_dir = resolve_gold_source_dir(source_dir, params)
    _, gold_state = load_state_from_dir(gold_dir)
    mean_dir = mean_evaluation_dir(output_dir)
    mean_source_dir = isnothing(seed_mean_dir) ? mean_dir : seed_mean_dir
    if !isnothing(seed_mean_dir) && !isdir(seed_mean_dir)
        println("Error: seed mean directory not found: $seed_mean_dir")
        return
    end
    mean_state, mean_count = mean_state_from_dir(mean_source_dir)
    if !isnothing(mean_state)
        state = state_with_mean_manifold(state, mean_state)
    end

    best_state, best_result = sparse_optimize(dataset, state, lambda; mean_state = mean_state, mean_run_count = mean_count, strict = strict)
    updated_mean_state, updated_mean_count = averaged_mean_state(mean_state, mean_count, best_state)
    updated_mean_cc = strict ? strict_mean_corr(best_result.manifold_raw, dataset.t, updated_mean_state) : bar_chart_corr(best_state, updated_mean_state)
    updated_strict_mean_cc = strict_mean_corr(best_result.manifold_raw, dataset.t, updated_mean_state)
    best_result = EvalResult(
        best_result.manifold_raw,
        best_result.manifold_fit,
        best_result.annual,
        best_result.final_fit,
        best_result.beta,
        best_result.cal_intercept,
        best_result.cal_slope,
        best_result.train_cc,
        best_result.validation_cc,
        best_result.test_cc,
        updated_mean_cc,
        updated_strict_mean_cc,
        best_result.mean_penalty_weight,
        best_result.mean_deviation,
        best_result.mean_penalty,
        best_result.complexity,
    )
    save_outputs(source_dir, output_dir, params, best_state, best_result, dataset; mean_dir = mean_dir, mean_run_count = updated_mean_count, strict = strict)
    save_mean_outputs(mean_dir, params, gold_dir, updated_mean_state, updated_mean_count)
    pruned_lpap, pruned_ks, pruned_params = pruned_summary(gold_state, best_state)

    println(
        "CC train=$(round(best_result.train_cc, digits=6)) " *
        "validation=$(round(best_result.validation_cc, digits=6)) " *
        "test=$(round(best_result.test_cc, digits=6)) " *
        "complexity=$(best_result.complexity)"
    )
    println("Observed months used=$(count(identity, dataset.observed_mask)) skipped_missing=$(length(dataset.observed_mask) - count(identity, dataset.observed_mask))")
    if !isnothing(seed_mean_dir)
        println("Seed mean loaded from $seed_mean_dir")
    end
    mean_label = strict ? "Strict mean manifold CC" : "Mean LPAP bar-chart CC"
    println("$mean_label=$(round(best_result.mean_manifold_cc, digits=6)) against $mean_dir (runs=$(updated_mean_count))")
    if !strict
        println("Strict mean manifold CC=$(round(best_result.strict_mean_manifold_cc, digits=6)) (--strict reference)")
    end
    println("Mean-manifold penalty: weight=$(round(best_result.mean_penalty_weight, digits=6)) deviation=$(round(best_result.mean_deviation, digits=6)) penalty=$(round(best_result.mean_penalty, digits=6))")
    println("Pruned vs gold ($gold_dir): lpap=[$pruned_lpap] active_ks=[$pruned_ks] params=[$pruned_params]")
    println("active_ks = retained modulation periods from the gold ltep/harm basis: [" * join(fmt_num.(active_ks(best_state)), ", ") * "]")
    println("beta = [intercept, manifold term, then sin/cos coefficients for active_ks in order]: " * beta_summary(best_result, best_state))
    println("Sparse evaluation saved to $output_dir")
end

main()
