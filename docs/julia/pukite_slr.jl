using JSON
using DelimitedFiles
using LinearAlgebra
using Statistics

# Pukite MSL Tidal Fit Algorithm
# Ported from pukite-slr.html to Julia
#
# Usage: julia pukite_slr.jl <experiment_dir> [lambda]
#
# This script reads from the following files in the experiment directory:
# - lt.exe.p: JSON parameters (bg, delA, delB, asym, year, mp, init, lpap, ltep, harm, etc.)
# - lt.exe.resp: Response file containing YEAR adjustment
# - lte_results.csv: Data file (date, source-model, observed, manifold-target, etc.)

const DEFAULT_IDATE = 1920.9
const HOLDOUT_START = 1940.0
const HOLDOUT_END = 1970.0

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
        if pd == 0
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
        # ampTime mirrors Ada Impulse_Amplify(Offset => bg, Ramp => year)
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
    # Forward integration
    for i in (si + 1):N
        ramp = m[i-1] >= 0 ? lag_c : -lag_c
        m[i] = forcing[i] + m[i-1] - ramp
    end
    # Backward integration
    for i in si:-1:2
        ramp = m[i] >= 0 ? lag_c : -lag_c
        m[i-1] = -forcing[i-1] + m[i] + ramp
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
    sxx = sum(x_vals.^2)
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

function corrcoef(x_vals, y_vals)
    return cor(collect(x_vals), collect(y_vals))
end

function parse_resp_year(resp_file)
    year_adj = 0.005 # Default
    if isfile(resp_file)
        lines = readlines(resp_file)
        for line in lines
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

function main()
    if length(ARGS) < 1
        println("Usage: julia pukite_slr.jl <experiment_dir> [lambda]")
        return
    end
    
    exp_dir = ARGS[1]
    lambda = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0
    
    p_file = joinpath(exp_dir, "lt.exe.p")
    resp_file = joinpath(exp_dir, "lt.exe.resp")
    csv_file = joinpath(exp_dir, "lte_results.csv")
    
    if !isfile(p_file) || !isfile(csv_file)
        println("Error: Required files not found in $exp_dir")
        return
    end
    
    params = JSON.parsefile(p_file)
    year_adj = parse_resp_year(resp_file)
    
    # Load data from CSV (no header)
    data = readdlm(csv_file, ',')
    t = data[:, 1]
    observed = data[:, 3]
    manifold_target = data[:, 4]
    
    # Extract parameters from JSON
    year_len = 365.2422484 + year_adj
    bg = params["bg"]
    del_a = params["delA"]
    del_b = params["delB"]
    asym = params["asym"]
    ann1 = params["ann1"]
    ann2 = params["ann2"]
    sem1 = params["sem1"]
    sem2 = params["sem2"]
    year_drift = params["year"]
    mp = params["mp"]
    init_val = params["init"]
    
    lpap_raw = params["lpap"]
    lpap_mat = zeros(Float64, length(lpap_raw), 3)
    for i in 1:length(lpap_raw)
        lpap_mat[i, :] = [Float64(x) for x in lpap_raw[i]]
    end
    
    ltep = params["ltep"]
    harm = params["harm"]
    ks = [ltep[1], ltep[2]]
    if length(ltep) > 1
        for h in harm
            push!(ks, ltep[2] * h)
        end
    end
    
    # Algorithm execution
    tide = build_tide(t, lpap_mat, year_len)
    forcing = build_forcing(t, tide, year_drift, bg, del_a, del_b, asym)
    manifold_raw = build_manifold(t, forcing, DEFAULT_IDATE, init_val, mp)
    annual = build_annual(t, ann1, ann2, sem1, sem2)
    
    # Affine calibration of manifold to match source target
    cal = solve_linear_map(manifold_raw, manifold_target)
    manifold_fit = apply_linear_map(manifold_raw, cal)
    
    # Multiple Linear Regression trained on the surrounding years, with
    # 1940-1970 treated as the hold-out validation window.
    y_resid = observed .- annual
    train_mask = (t .< HOLDOUT_START) .| (t .> HOLDOUT_END)
    holdout_mask = .!train_mask
    beta = solve_mlr(manifold_fit[train_mask], y_resid[train_mask], ks, lambda)
    fit_resid = apply_beta(manifold_fit, beta, ks)
    final_fit = fit_resid .+ annual
    train_cc = corrcoef(observed[train_mask], final_fit[train_mask])
    holdout_cc = corrcoef(observed[holdout_mask], final_fit[holdout_mask])
    
    # Save results to CSV in the experiment directory
    output_file = joinpath(exp_dir, "julia_results.csv")
    # Columns: Date, Raw Manifold, Calibrated Manifold, Final Fit
    writedlm(output_file, hcat(t, manifold_raw, manifold_fit, final_fit), ',')
    println("CC train=$(round(train_cc, digits=6)) holdout=$(round(holdout_cc, digits=6))")
    println("Algorithm completed. Results saved to $output_file")
end

main()
