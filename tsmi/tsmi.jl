#!/usr/bin/env julia

using DelimitedFiles
using JSON3
using Dates
using LinearAlgebra
using Statistics
using Optim
using DSP

# ---------- I/O helpers ----------

function read_timeseries(path::String)
    data = readdlm(path)  # whitespace-delimited
    t = data[:, 1]
    I = data[:, 2]
    return t, I
end

function read_params(path::String)
    txt = read(path, String)
    j = JSON3.read(txt, Dict)
    tides = j["tides"]          # array of objects with period, amplitude, phase
    damping = j["damping"]      # object with zeta, omega0
    kappa = collect(-2000:2000) # "Large faster LTE modulation"
    annual_phase = Float64(get(j, "annual_phase", 0.0))
    
    # Read initial condition parameters (default to 0.0)
    ic = get(j, "initial_condition", Dict("A"=>0.0, "B"=>0.0))
    ic_A = Float64(get(ic, "A", 0.0))
    ic_B = Float64(get(ic, "B", 0.0))
    
    # Read ramp slope if available, else default to 0.001 (small decay per step)
    ramp_slope = Float64(get(j, "ramp_slope", 0.001))
    
    # Read DC offset
    dc_offset = Float64(get(j, "dc_offset", 0.0))
    
    # Read Impulse Sigma (default to 1.0/78.0 approx 0.0128)
    impulse_sigma = Float64(get(j, "impulse_sigma", 1.0/78.0))
    
    # Read Primary Impulse Amplitude (default 1.0)
    primary_amp = Float64(get(j, "primary_impulse_amp", 1.0))
    
    # Read IC Reference Index (default 1)
    ic_ref_index = Int(get(j, "ic_ref_index", 1))
    
    # Read Drift Pivot Time (optional, default nothing)
    t_pivot = nothing
    if haskey(j, "drift_pivot_time")
        t_pivot = Float64(j["drift_pivot_time"])
    end
    
    return tides, damping, kappa, annual_phase, ic_A, ic_B, ramp_slope, dc_offset, impulse_sigma, primary_amp, ic_ref_index, t_pivot, j
end

function archive_params(path::String)
    ts = Dates.format(now(), "yyyymmdd_HHMMSS")
    backup = "$(path).bak_$ts"
    cp(path, backup; force=true)
    return backup
end

function write_params(path::String, json_obj)
    open(path, "w") do io
        JSON3.write(io, json_obj)
    end
    # Use python to pretty print since JSON3 pretty print is limited
    run(`python -c "import json; d=json.load(open('$path')); open('$path', 'w').write(json.dumps(d, indent=4))"`)
end

function sparsify_coeffs(kappa, c_dense)
    sparse_list = []
    for (i, val) in enumerate(c_dense)
        if abs(val) > 1e-9
            push!(sparse_list, Dict("k" => kappa[i], "re" => real(val), "im" => imag(val)))
        end
    end
    return sparse_list
end

function load_active_coeffs(json_obj)
    # Returns Vector{Tuple{Int, ComplexF64}}
    active_coeffs = Tuple{Int, ComplexF64}[]
    
    if !haskey(json_obj, "coefficients"); return active_coeffs; end
    
    data = json_obj["coefficients"]
    
    if data isa Vector
        # New sparse format: list of objects
        for item in data
            # Handle JSON3 object access
            k = Int(item["k"])
            val = ComplexF64(Float64(item["re"]), Float64(item["im"]))
            push!(active_coeffs, (k, val))
        end
    elseif data isa Dict
        # Old dense format: dict with arrays
        if haskey(data, "kappa")
            kappa_stored = Int.(data["kappa"])
            re = Float64.(data["real"])
            im = Float64.(data["imag"])
            for (i, k) in enumerate(kappa_stored)
                val = ComplexF64(re[i], im[i])
                if abs(val) > 1e-9
                    push!(active_coeffs, (k, val))
                end
            end
        end
    end
    return active_coeffs
end

# ---------- Physics: manifold construction ----------

function build_tidal_signal(t, tides, year_len)
    # t is in years. Convert to days using the variable year length.
    t_days = t .* year_len
    sig = zeros(eltype(t), length(t))
    for td in tides
        P  = Float64(td["period"])
        A  = Float64(td["amplitude"])
        φ  = Float64(td["phase"])
        ω  = 2π / P
        @. sig += A * sin(ω * t_days + φ)
    end
    return sig
end

function build_annual_comb(t, primary_amp, phi, sec_amp, year_len, drift_linear, drift_quad, impulse_sigma; t_pivot=nothing)
    # Gaussian comb with period 365.0 days + drifts
    # t is in years (approx 365.25 days)
    # P_comb = 365.0 / 365.25 approx 0.999... in 't' units?
    # No, let's keep it simple: 
    # The 'nominal' year length in the model is year_len.
    # The impulse happens every 365.0 days.
    # So P_comb = 365.0 / year_len (in units of 't')
    
    P_comb = 365.0 / year_len
    
    # Sigma for Gaussian Impulse
    # User feedback indicates "plateaus" should be smooth. 
    # A wide impulse (1/12) causes tidal wiggles to leak into the plateau.
    # Narrowing to ~1 week (1/52) to create cleaner steps/plateaus.
    σ = impulse_sigma # default 1.0/78.0
    
    comb = zeros(eltype(t), length(t))
    
    # Range of k to cover t
    # Center time for quadratic drift? 
    # If t_pivot is provided, use it. Otherwise use mean(t).
    # t_pivot is critical for consistent application of drift parameters across different time intervals.
    t_mean = isnothing(t_pivot) ? mean(t) : t_pivot
    
    k_min = floor(Int, minimum(t) / P_comb) - 2
    k_max = ceil(Int, maximum(t) / P_comb) + 2
    
    for k in k_min:k_max
        # Nominal time
        t_nom = k * P_comb + phi
        
        # Add drifts
        # Linear drift is already covered by P_comb vs year_len mismatch?
        # User says "slight linear drift AND quadratic drift in impulse time"
        # So: t_actual = t_nom + L * (t_nom - t_mean) + Q * (t_nom - t_mean)^2
        
        dt = t_nom - t_mean
        t_actual = t_nom + drift_linear * dt + drift_quad * dt^2
        
        # Primary Impulse (Variable Amplitude)
        @. comb += primary_amp * exp( - (t - t_actual)^2 / (2 * σ^2) )
        
        # Secondary Impulse
        # Assume secondary impulse shifts with the primary
        t_sec = t_actual + 0.5 * P_comb 
        @. comb += sec_amp * exp( - (t - t_sec)^2 / (2 * σ^2) )
    end
    return comb
end

function build_manifold(t, tides, ramp_slope, primary_amp, phi, sec_amp, year_len, ic_A, ic_B, drift_linear, drift_quad, 
                        comp_ann_amp, comp_ann_phase, comp_semi_amp, comp_semi_phase, dc_offset, impulse_sigma; include_additive=true, ref_index=1, t_pivot=nothing)
    # Spin-up: Run for 20 years prior to t[1] to eliminate transient
    
    t_start = t[1]
    dt = t[2] - t[1] 
    if dt <= 0; dt = 1.0/365.242; end
    
    t_spinup = range(t_start - 40.0, stop=t_start - dt, step=dt) 
    
    # Base Multiplicative (Impulse * Tides)
    # Note: For spin-up, we should ideally use the same t_pivot, but spin-up is just for settling
    # Use the passed t_pivot if available for consistency
    mult_spinup = build_annual_comb(t_spinup, primary_amp, phi, sec_amp, year_len, drift_linear, drift_quad, impulse_sigma; t_pivot=t_pivot) .* build_tidal_signal(t_spinup, tides, year_len)
    
    # Additive Compensatory Sinusoids
    omega_ann = 2π
    if include_additive
        comp_spinup = @. comp_ann_amp * sin(omega_ann * t_spinup + comp_ann_phase) + 
                         comp_semi_amp * sin(2 * omega_ann * t_spinup + comp_semi_phase)
        base_spinup = mult_spinup .+ comp_spinup
    else
        base_spinup = mult_spinup
    end
    
    # Spin-up integration (still useful to establish trend if needed, but we override final value)
    val = 0.0 
    for i in eachindex(t_spinup)
        ramp = copysign(ramp_slope, val)
        val = base_spinup[i] + val - ramp
    end
    
    # Main run
    mult = build_annual_comb(t, primary_amp, phi, sec_amp, year_len, drift_linear, drift_quad, impulse_sigma; t_pivot=t_pivot) .* build_tidal_signal(t, tides, year_len)
    
    if include_additive
        comp = @. comp_ann_amp * sin(omega_ann * t + comp_ann_phase) + 
                  comp_semi_amp * sin(2 * omega_ann * t + comp_semi_phase)
        base = mult .+ comp
    else
        base = mult
    end
    
    M_forced = zeros(eltype(t), length(t))
    
    # ---------------------------------------------------------
    # Integration Logic (Forward/Backward from ref_index)
    # ---------------------------------------------------------
    
    if ref_index < 1 || ref_index > length(t)
        ref_index = 1
    end
    
    # Forward Integration (from ref_index to end)
    # ic_A is defined as the accumulator value BEFORE step ref_index.
    # I.e., X[ref_index-1] = ic_A
    
    current_val = ic_A
    
    for i in ref_index:length(t)
        ramp = copysign(ramp_slope, current_val)
        current_val = base[i] + current_val - ramp
        M_forced[i] = current_val + dc_offset
    end
    
    # Backward Integration (from ref_index-1 down to 1)
    if ref_index > 1
        # Set the value at ref_index-1 explicitly (it is our starting condition)
        M_forced[ref_index-1] = ic_A + dc_offset
        
        current_val_back = ic_A
        
        # We need to compute X[i-1] from X[i]
        # Step i: X[i] = X[i-1] + base[i] - ramp(X[i-1])
        # We iterate i descending.
        # We know X[i] (current_val_back). We want X[i-1] (prev_val).
        # But wait, M indices match t indices.
        # We have X[ref_index-1]. We want X[ref_index-2].
        # Step ref_index-1 produced X[ref_index-1] from X[ref_index-2].
        # So we invert step k = ref_index-1, ref_index-2, ... 2.
        
        for i in (ref_index-1):-1:2
            # Inverting step i
            # X[i] = X[i-1] + base[i] - ramp(X[i-1])
            # Z = X[i] - base[i] = X[i-1] - ramp(X[i-1])
            
            Z = current_val_back - base[i]
            
            # Apply inverse ramp correction
            # Using Z as proxy for X[i-1] to determine sign
            ramp = copysign(ramp_slope, Z)
            prev_val = Z + ramp
            
            M_forced[i-1] = prev_val + dc_offset
            current_val_back = prev_val
        end
    end
    
    return M_forced
end



# ---------- Sparsity Helpers ----------

# ---------- OMP Core ----------

function run_omp(M, I, kappa, t; max_atoms=18)
    # Construct Manifold Dictionary
    n_kappa = length(kappa)
    n_add = 4 # 2 annual + 2 semi-annual
    Phi = Matrix{ComplexF64}(undef, length(I), n_kappa + n_add)
    
    # Manifold terms
    @inbounds for (j, k_val) in enumerate(kappa)
        @. Phi[:, j] = exp(1im * k_val * M)
    end
    
    # Additive terms (Annual + Semi-Annual)
    # Indices: n_kappa + 1..4
    # Period 1.0 (Annual) -> 2pi * t
    omega_ann = 2π
    @. Phi[:, n_kappa + 1] = cos(omega_ann * t)
    @. Phi[:, n_kappa + 2] = sin(omega_ann * t)
    @. Phi[:, n_kappa + 3] = cos(2 * omega_ann * t)
    @. Phi[:, n_kappa + 4] = sin(2 * omega_ann * t)
    
    residual = copy(ComplexF64.(I))
    selected_indices = Int[]
    c_full = zeros(ComplexF64, n_kappa + n_add)
    
    best_corr = 0.0
    
    for iter in 1:max_atoms
        projections = Phi' * residual
        mags = abs.(projections)
        
        if !isempty(selected_indices)
            mags[selected_indices] .= -1.0
        end
        
        best_idx = argmax(mags)
        if mags[best_idx] <= 1e-9
            break
        end
        
        push!(selected_indices, best_idx)
        
        Phi_sub = Phi[:, selected_indices]
        # Ridge for stability
        lambda = 1e-6
        H = Phi_sub' * Phi_sub
        for i in 1:size(H,1); H[i,i] += lambda; end
        c_sub = H \ (Phi_sub' * ComplexF64.(I))
        
        Im_sub = Phi_sub * c_sub
        residual = ComplexF64.(I) - Im_sub
        
        Im_real = real.(Im_sub)
        r = cor(I, Im_real)
        if r > best_corr
            best_corr = r
        end
        
        c_full[:] .= 0.0
        c_full[selected_indices] = c_sub
    end
    
    return best_corr, c_full
end

# ---------- Optimization ----------                  0.05

function mutate_relative(val::Float64; scale::Float64=0.01, neg_prob::Float64=0.1, floor_val::Float64=1e-9)
    factor = 1.0 + scale * randn()
    if rand() < neg_prob
        factor = -factor
    end
    result = val * factor
    if abs(result) < floor_val
        result = copysign(floor_val, rand() < 0.5 ? 1.0 : -1.0)
    end
    return result
end

function optimize_manifold_random(t, I, json_obj, ic_A, ic_B, kappa; manifold_only=false, stop_threshold=nothing, ref_index=1, t_pivot=nothing)
    println("Starting Canonical Manifold Optimization (Random Search - Linear/Quad Drift)...")
    if manifold_only
        println("  (Mode: Manifold-Only Fit)")
    end
    if stop_threshold !== nothing
        println("  (Stop Threshold: Correlation >= $stop_threshold)")
    end
    println("  (Integration Reference Index: $ref_index)")
    if t_pivot !== nothing
        println("  (Drift Pivot Time: $t_pivot)")
    end
    
    # Optimization Target is always the input (t, I)
    I_target = I
    t_target = t

    # Current Params
    tides = json_obj["tides"]
    current_tides = deepcopy(tides)
    current_ramp = Float64(get(json_obj, "ramp_slope", 0.001))
    current_phi = Float64(get(json_obj, "annual_phase", 0.0))
    current_year_len = Float64(get(json_obj, "year_length", 365.242))
    current_primary_amp = Float64(get(json_obj, "primary_impulse_amp", 1.0))
    current_sec_amp = Float64(get(json_obj, "secondary_impulse_amp", 0.0))
    current_d_lin = Float64(get(json_obj, "drift_linear", 0.0))
    current_d_quad = Float64(get(json_obj, "drift_quad", 0.0))
    current_ca_amp = Float64(get(json_obj, "comp_annual_amp", 0.0))
    current_ca_phase = Float64(get(json_obj, "comp_annual_phase", 0.0))
    current_csa_amp = Float64(get(json_obj, "comp_semi_amp", 0.0))
    current_csa_phase = Float64(get(json_obj, "comp_semi_phase", 0.0))
    current_dc = Float64(get(json_obj, "dc_offset", 0.0))
    current_sigma = Float64(get(json_obj, "impulse_sigma", 1.0/78.0))
    
    # Use input initial conditions as starting point
    current_ic_A = ic_A
    current_ic_B = ic_B
    
    best_tides = deepcopy(current_tides)
    best_ramp = current_ramp
    best_phi = current_phi
    best_year_len = current_year_len
    best_primary_amp = current_primary_amp
    best_sec_amp = current_sec_amp
    best_d_lin = current_d_lin
    best_d_quad = current_d_quad
    best_ca_amp = current_ca_amp
    best_ca_phase = current_ca_phase
    best_csa_amp = current_csa_amp
    best_csa_phase = current_csa_phase
    best_dc = current_dc
    best_sigma = current_sigma
    best_ic_A = current_ic_A
    best_ic_B = current_ic_B
    
    # Objective function wrapper
    sigma_I = std(I_target) # Compute once outside
    function calc_objective(M_trial, I_ref)
        if manifold_only
            # For manifold fitting, we want to match the absolute amplitude and shape.
            # Avoid linear scaling (slope/intercept) to force amplitude convergence.
            
            # Calculate Absolute RMSE
            diff = M_trial .- I_ref
            mse = mean(diff.^2)
            rmse_abs = sqrt(mse)
            
            # Calculate Correlation
            r = cor(M_trial, I_ref)
            
            # Combined Metric:
            # We want high correlation (r -> 1) and low RMSE (rmse_abs -> 0).
            # Normalize RMSE by the standard deviation of the target signal to make it comparable.
            nrmse = rmse_abs / (sigma_I + 1e-9)
            
            # Score formulation:
            # If amplitude is off by factor of 2 (e.g. 0.03 vs 0.06), RMSE will be large.
            # Let's weight RMSE heavily to force amplitude matching.
            # Score = Correlation - Weight * NRMSE
            # If Weight is too high, it might sacrifice correlation for DC offset?
            # But M and I should be aligned.
            
            score = r - 0.1 * nrmse  # 2
            return score
        else
            # Correlation with OMP reconstruction (General Case)
            # If I_ref is forcing, only need 1 atom to correlate well.
            # But generally we want to fit I_ref.
            # Use simpler OMP for speed? 
            r, _ = run_omp(M_trial, I_ref, kappa, t_target; max_atoms=18)
            return r
        end
    end

    # Determine if we include additive terms in manifold construction
    # If manifold_only is true, we EXCLUDE additive terms (ramp only) per requirement
    inc_add = !manifold_only

    # Initial Baseline
    M = build_manifold(t_target, current_tides, current_ramp, current_primary_amp, current_phi, current_sec_amp, current_year_len, current_ic_A, current_ic_B, current_d_lin, current_d_quad,
                       current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase, current_dc, current_sigma; include_additive=inc_add, ref_index=ref_index, t_pivot=t_pivot)
    r_best = calc_objective(M, I_target)
    
    println("Initial Correlation: $r_best")

    ################################################
    
    # Pre-Optimization: Grid Search for Initial Conditions (ic_A)
    # The initial condition of the integrator affects the transient at the start of the time-series.
    # Scan a reasonable range (e.g., -2.0 to 2.0) to match the starting level of the target data.
    println("Running Grid Search for Initial Condition A...")
    best_ic_scan_A = current_ic_A
    r_scan_best_ic = r_best
    
    # Range based on typical normalized anomaly data (-2.0 to 2.0)
    ic_scan_steps = range(-4.0, stop=4.0, length=200) 
    for ic_val in ic_scan_steps
        M_scan = build_manifold(t_target, current_tides, current_ramp, current_primary_amp, current_phi, current_sec_amp, current_year_len, ic_val, current_ic_B, current_d_lin, current_d_quad,
                               current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase, current_dc, current_sigma; include_additive=inc_add, ref_index=ref_index, t_pivot=t_pivot)
        r_scan = calc_objective(M_scan, I_target)
        if r_scan > r_scan_best_ic
            r_scan_best_ic = r_scan
            best_ic_scan_A = ic_val
        end
    end
    
    if r_scan_best_ic > r_best
        println("  Found better Initial Condition A: $best_ic_scan_A (r=$r_scan_best_ic vs previous $r_best)")
        # current_ic_A = best_ic_scan_A
        # best_ic_A = best_ic_scan_A
        # r_best = r_scan_best_ic
    else
        println("  Initial Condition A retained.")
    end

    ################################

    # Pre-Optimization: Grid Search for Annual Phase
    # The annual phase (impulse timing) can have local optima due to narrow impulse width.
    # Scan through one full period (~1.0 year) to find the best basin.
    if ref_index == 1
        println("Running Grid Search for Annual Phase...")
    best_phase_scan = current_phi
    r_scan_best = r_best
    
    # 24 steps (twice per month)
    scan_steps = range(0.0, stop=1.0, length=25)[1:end-1] 
    for p_val in scan_steps
        M_scan = build_manifold(t_target, current_tides, current_ramp, current_primary_amp, p_val, current_sec_amp, current_year_len, current_ic_A, current_ic_B, current_d_lin, current_d_quad,
                               current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase, current_dc, current_sigma; include_additive=inc_add, ref_index=ref_index, t_pivot=t_pivot)
        r_scan = calc_objective(M_scan, I_target)
        if r_scan > r_scan_best
            r_scan_best = r_scan
            best_phase_scan = p_val
        end
    end
    
    if r_scan_best > r_best
        println("  Found better phase start: $best_phase_scan (r=$r_scan_best vs initial $r_best)")
        current_phi = best_phase_scan
        best_phi = best_phase_scan
        r_best = r_scan_best
    else
        println("  Initial phase retained.")
    end

    # Pre-Optimization: Grid Search for Initial Conditions (ic_A)
    println("Running Grid Search for Initial Condition A...")
    best_ic_scan_A = current_ic_A
    r_scan_best_ic = r_best
    
    # Range based on typical normalized anomaly data (-2.0 to 2.0)
    ic_scan_steps = range(-2.0, stop=2.0, length=41) 
    for ic_val in ic_scan_steps
        M_scan = build_manifold(t_target, current_tides, current_ramp, current_primary_amp, current_phi, current_sec_amp, current_year_len, ic_val, current_ic_B, current_d_lin, current_d_quad,
                               current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase, current_dc, current_sigma; include_additive=inc_add, ref_index=ref_index, t_pivot=t_pivot)
        r_scan = calc_objective(M_scan, I_target)
        if r_scan > r_scan_best_ic
            r_scan_best_ic = r_scan
            best_ic_scan_A = ic_val
        end
    end
    
    if r_scan_best_ic > r_best
        println("  Found better Initial Condition A: $best_ic_scan_A (r=$r_scan_best_ic vs previous $r_best)")
        current_ic_A = best_ic_scan_A
        best_ic_A = best_ic_scan_A
        r_best = r_scan_best_ic
    else
        println("  Initial Condition A retained.")
    end
    else
        println("Skipping Grid Searches (Reference Index != 1)")
    end

    ################################################


    # Loop
    # Increase iterations for manifold-only mode to find global optimum
    n_iter = manifold_only ? 250000 : 5000
    
    # Adaptive Temperature Schedule
    # If starting from a good fit (high r_best), start colder to avoid destroying progress.
    # Also ensures we don't freeze too early by setting alpha based on n_iter.
    err = max(0.001, 1.0 - r_best)
    T_start = min(0.1, err * 0.5) 
    T_min = 1e-7
    
    # Calculate alpha to decay from T_start to T_min over n_iter
    alpha = (T_min / T_start)^(1.0 / n_iter)
    
    T = T_start
    println("  Adaptive Schedule: T_start=$(round(T, digits=6)), alpha=$(round(alpha, digits=8)) for $n_iter iters")
    
    idx_27 = findfirst(x -> abs(Float64(x["period"]) - 27.2122) < 0.1, current_tides)
    
    # Current state
    curr_tides = deepcopy(best_tides)
    curr_ramp = best_ramp
    curr_phi = best_phi # Ensure curr_phi starts at the best found phase
    curr_year_len = best_year_len
    curr_primary_amp = best_primary_amp
    curr_sec_amp = best_sec_amp
    curr_d_lin = best_d_lin
    curr_d_quad = best_d_quad
    curr_ca_amp = best_ca_amp
    curr_ca_phase = best_ca_phase
    curr_csa_amp = best_csa_amp
    curr_csa_phase = best_csa_phase
    curr_dc = best_dc
    curr_ic_A = best_ic_A
    curr_ic_B = best_ic_B
    curr_sigma = best_sigma
    curr_r = r_best
    
    for i in 1:n_iter
        cand_tides = deepcopy(curr_tides)
        cand_ramp = curr_ramp
        cand_phi = curr_phi
        cand_year_len = curr_year_len
        cand_primary_amp = curr_primary_amp
        cand_sec_amp = curr_sec_amp
        cand_d_lin = curr_d_lin
        cand_d_quad = curr_d_quad
        cand_ca_amp = curr_ca_amp
        cand_ca_phase = curr_ca_phase
        cand_csa_amp = curr_csa_amp
        cand_csa_phase = curr_csa_phase
        cand_dc = curr_dc
        cand_ic_A = curr_ic_A
        cand_ic_B = curr_ic_B
        cand_sigma = curr_sigma
        
        # Mutate Tides
        for (idx, td) in enumerate(cand_tides)
            if rand() < 0.05 || (idx == idx_27 && rand() < 0.2)
                td["amplitude"] = mutate_relative(Float64(td["amplitude"]))
                td["phase"] = Float64(td["phase"]) + randn() * 0.1
            end
        end
        
        # Mutate Global Params
        if rand() < 0.2; cand_ramp = mutate_relative(cand_ramp; scale=0.1, neg_prob=0.0); end
        if rand() < 0.2; cand_primary_amp = mutate_relative(cand_primary_amp); end
        if rand() < 0.2; cand_phi += randn() * 0.05; end
        if rand() < 0.2; cand_year_len = mutate_relative(cand_year_len; scale=0.001, neg_prob=0.0); end
        if rand() < 0.2; cand_sec_amp = mutate_relative(cand_sec_amp); end
        
        # Mutate Drifts (Very small steps)
        if rand() < 0.3; cand_d_lin = mutate_relative(cand_d_lin); end
        if rand() < 0.3; cand_d_quad = mutate_relative(cand_d_quad); end
        
        # Mutate Compensatory Terms (Only if included)
        if inc_add
            if rand() < 0.2
                cand_ca_amp = mutate_relative(cand_ca_amp)
                cand_ca_phase += randn() * 0.1
            end
            if rand() < 0.2
                cand_csa_amp = mutate_relative(cand_csa_amp)
                cand_csa_phase += randn() * 0.1
            end
        end

        # Mutate Initial Conditions and DC
        if rand() < 0.2
            cand_ic_A = mutate_relative(cand_ic_A)
            cand_ic_B = mutate_relative(cand_ic_B)
        end
        if rand() < 0.2
            # DC offset mutation (additive, relative to signal scale)
            cand_dc += randn() * (sigma_I * 0.005)
        end
        
        # Mutate Impulse Width (Sigma)
        if rand() < 0.2
            cand_sigma = mutate_relative(cand_sigma)
        end
        
        # Evaluate
        M_cand = build_manifold(t_target, cand_tides, cand_ramp, cand_primary_amp, cand_phi, cand_sec_amp, cand_year_len, cand_ic_A, cand_ic_B, cand_d_lin, cand_d_quad,
                                cand_ca_amp, cand_ca_phase, cand_csa_amp, cand_csa_phase, cand_dc, cand_sigma; include_additive=inc_add, ref_index=ref_index, t_pivot=t_pivot)
        r_cand = calc_objective(M_cand, I_target)
        
        # Metropolis-Hastings
        delta = r_cand - curr_r
        if delta > 0 || rand() < exp(delta * 10.0 / T) 
            curr_tides = cand_tides
            curr_ramp = cand_ramp
            curr_primary_amp = cand_primary_amp
            curr_phi = cand_phi
            curr_year_len = cand_year_len
            curr_sec_amp = cand_sec_amp
            curr_d_lin = cand_d_lin
            curr_d_quad = cand_d_quad
            curr_ca_amp = cand_ca_amp
            curr_ca_phase = cand_ca_phase
            curr_csa_amp = cand_csa_amp
            curr_csa_phase = cand_csa_phase
            curr_dc = cand_dc
            curr_sigma = cand_sigma
            curr_ic_A = cand_ic_A
            curr_ic_B = cand_ic_B
            curr_r = r_cand
            
            if r_cand > r_best
                println("Iter $i: New Best r = $r_cand (was $r_best)")
                r_best = r_cand
                best_tides = deepcopy(cand_tides)
                best_ramp = cand_ramp
                best_primary_amp = cand_primary_amp
                best_phi = cand_phi
                best_year_len = cand_year_len
                best_sec_amp = cand_sec_amp
                best_d_lin = cand_d_lin
                best_d_quad = cand_d_quad
                best_ca_amp = cand_ca_amp
                best_ca_phase = cand_ca_phase
                best_csa_amp = cand_csa_amp
                best_csa_phase = cand_csa_phase
                best_dc = cand_dc
                best_sigma = cand_sigma
                best_ic_A = cand_ic_A
                best_ic_B = cand_ic_B
                
                if stop_threshold !== nothing && r_best >= stop_threshold
                    println("  Hit target correlation $stop_threshold at iter $i. Stopping early.")
                    break
                end
            end
        end
        
        T *= alpha
        if i % 500 == 0
            println("  Iter $i... Temp: $(round(T, digits=4)), Curr: $(round(curr_r, digits=4)), Best: $(round(r_best, digits=4))")
        end
    end
    
    return best_tides, best_ramp, best_primary_amp, best_phi, best_sec_amp, best_year_len, best_d_lin, best_d_quad, best_ca_amp, best_ca_phase, best_csa_amp, best_csa_phase, best_dc, best_sigma, best_ic_A, best_ic_B, r_best
end

# ---------- Main ----------

function regression_test_integration(t, tides, ramp_slope, primary_amp, phi, sec_amp, year_len, drift_linear, drift_quad, impulse_sigma, dc_offset, json_obj)
    println("\nRunning Regression Test for Arbitrary Start Point...")
    
    # 1. Run Standard Forward (ref_index = 1)
    # We need a baseline ic_A. Let's use 0.0 or a typical value.
    ic_A_base = 0.5
    ic_B_base = 0.0 # Unused
    
    # Ensure additive is consistent
    ca_amp = 0.0; ca_ph = 0.0; csa_amp = 0.0; csa_ph = 0.0
    
    M1 = build_manifold(t, tides, ramp_slope, primary_amp, phi, sec_amp, year_len, ic_A_base, ic_B_base, drift_linear, drift_quad,
                       ca_amp, ca_ph, csa_amp, csa_ph, dc_offset, impulse_sigma; include_additive=false, ref_index=1)
                       
    # 2. Pick a mid-point
    mid_idx = length(t) ÷ 2
    
    # 3. Get value at mid-point (remove DC to get accumulator value)
    # ic_A is defined as the value BEFORE the ref_index step.
    # So if ref_index = mid_idx, we need the value at mid_idx-1.
    ic_mid = M1[mid_idx-1] - dc_offset
    
    # 4. Run with ref_index = mid_idx
    M2 = build_manifold(t, tides, ramp_slope, primary_amp, phi, sec_amp, year_len, ic_mid, ic_B_base, drift_linear, drift_quad,
                       ca_amp, ca_ph, csa_amp, csa_ph, dc_offset, impulse_sigma; include_additive=false, ref_index=mid_idx)
                       
    # 5. Compare
    diff = M1 .- M2
    max_diff = maximum(abs.(diff))
    
    # Check forward and backward parts separately
    diff_forward = maximum(abs.(diff[mid_idx:end]))
    diff_backward = maximum(abs.(diff[1:mid_idx-1]))
    
    println("  Max Difference Forward (Ref to End): $diff_forward")
    println("  Max Difference Backward (Start to Ref): $diff_backward")
    
    if max_diff < 1e-9
        println("  [PASS] Integration is fully reversible/consistent.")
    else
        if diff_forward < 1e-9
            println("  [PARTIAL] Forward integration consistent. Backward integration has divergence ($diff_backward).")
            println("            (Note: Backward divergence is expected due to dissipative ramp function)")
        else
            println("  [FAIL] Forward integration inconsistency detected ($diff_forward).")
        end
    end
    
    # 6. Verify Logic with Zero Slope (Should be perfect)
    println("  Verifying Indexing Logic (Zero Slope)...")
    M1_z = build_manifold(t, tides, 0.0, primary_amp, phi, sec_amp, year_len, ic_A_base, ic_B_base, drift_linear, drift_quad,
                         ca_amp, ca_ph, csa_amp, csa_ph, dc_offset, impulse_sigma; include_additive=false, ref_index=1)
    ic_mid_z = M1_z[mid_idx-1] - dc_offset
    M2_z = build_manifold(t, tides, 0.0, primary_amp, phi, sec_amp, year_len, ic_mid_z, ic_B_base, drift_linear, drift_quad,
                         ca_amp, ca_ph, csa_amp, csa_ph, dc_offset, impulse_sigma; include_additive=false, ref_index=mid_idx)
    diff_z = maximum(abs.(M1_z .- M2_z))
    if diff_z < 1e-9
        println("  [PASS] Zero-slope logic check passed (diff=$diff_z). Indices are correct.")
    else
        println("  [FAIL] Zero-slope logic check failed (diff=$diff_z). Indexing error persists.")
    end
end

function main(ts_path::String, json_path::String; skip_optim=false, manifold_only=false, final_mod_path::Union{String, Nothing}=nothing, staged_path::Union{String, Nothing}=nothing)
    t, I = read_timeseries(ts_path)
    tides, damping, kappa, annual_phase, ic_A, ic_B, ramp_slope, dc_param, impulse_sigma, primary_amp, ic_ref_index, t_pivot_in, json_obj = read_params(json_path)

    # ------------------------------------------------------------------
    # Automatic Initial Condition & Pivot Bookkeeping
    # ------------------------------------------------------------------
    
    # 1. Pivot Time
    # If t_pivot_in is missing, set it to the mean of current t (Training Phase logic)
    # If t_pivot_in is present, use it (Testing Phase logic).
    # But wait, if this IS the training phase (e.g. first run), we want to establish it.
    # If t_pivot_in is nothing, we establish it.
    if t_pivot_in === nothing
        t_pivot_in = mean(t)
        println("  [Bookkeeping] Establishing new Drift Pivot Time: $t_pivot_in")
        json_obj["drift_pivot_time"] = t_pivot_in
    else
        println("  [Bookkeeping] Using existing Drift Pivot Time: $t_pivot_in")
    end

    # 2. Initial Condition Reference Index
    # Check if we have a stored reference TIME.
    stored_ref_time = get(json_obj, "initial_condition_ref_time", nothing)
    
    if stored_ref_time !== nothing
        stored_ref_time = Float64(stored_ref_time)
        # We have a stored reference time. Check if current ic_ref_index matches it in current t.
        # If t[ic_ref_index] != stored_ref_time, we have a mismatch (new file?).
        
        current_time_at_idx = (ic_ref_index >= 1 && ic_ref_index <= length(t)) ? t[ic_ref_index] : nothing
        
        if current_time_at_idx === nothing || abs(current_time_at_idx - stored_ref_time) > 1e-4
            println("  [Bookkeeping] Detected time-series change or index mismatch.")
            println("    Stored Ref Time: $stored_ref_time")
            println("    Current Time at Index $ic_ref_index: $current_time_at_idx")
            
            # Find closest index in current t to stored_ref_time
            # Search t for closest value
            dists = abs.(t .- stored_ref_time)
            min_dist, new_idx = findmin(dists)
            
            if min_dist < 0.5 # Within half a year (reasonable for monthly/daily data)
                println("    Found matching time in new series at index $new_idx (diff $min_dist). Updating ref_index.")
                ic_ref_index = new_idx
                json_obj["ic_ref_index"] = new_idx
            else
                println("    [WARNING] Stored reference time $stored_ref_time is far from any value in current series (closest diff $min_dist).")
                # Fallback: Transport
                println("    Attempting to transport Initial Condition from $stored_ref_time to $(t[1])...")
                
                dt = 1.0/365.242
                if length(t) > 1; dt = t[2] - t[1]; end
                
                # Create a time vector covering the gap
                t_gap = collect(range(stored_ref_time, stop=t[1], step=dt))
                # Remove last point if it overlaps with t[1]
                if length(t_gap) > 0 && abs(t_gap[end] - t[1]) < 1e-6
                     pop!(t_gap) 
                end
                
                if length(t_gap) > 0
                    # Run physics over gap
                    drift_L = get(json_obj, "drift_linear", 0.0)
                    drift_Q = get(json_obj, "drift_quad", 0.0)
                    
                    # We need to build the manifold (integration)
                    # We treat ref_index=1 (at stored_ref_time)
                    # We integrate forward.
                    M_gap = build_manifold(t_gap, tides, ramp_slope, primary_amp, annual_phase, 0.0, 365.242, ic_A, ic_B, 
                                          drift_L, drift_Q, 0.0, 0.0, 0.0, 0.0, dc_param, impulse_sigma; 
                                          include_additive=false, ref_index=1, t_pivot=t_pivot_in)
                    
                    # The final value M_gap[end] + step -> t[1]
                    # The integrator state at the END of step N (M_gap[end]) is the value before step N+1.
                    # Since t_gap[end] is one step before t[1], M_gap[end] is the correct IC for t[1].
                    # Note: M_gap contains DC offset. Remove it.
                    final_state = M_gap[end] - dc_param
                    
                    println("    Transported IC_A from $ic_A (at $stored_ref_time) to $final_state (at $(t_gap[end])).")
                    
                    ic_A = final_state
                    ic_ref_index = 1
                    
                    # Update JSON object so we use the new IC
                    json_obj["initial_condition"]["A"] = ic_A
                    json_obj["ic_ref_index"] = 1
                    json_obj["initial_condition_ref_time"] = t[1]
                    println("    [Bookkeeping] Updated Initial Condition A to transported value.")
                end
            end
        end
    else
        # No stored ref time. First run?
        # Store the current time at ref_index
        if ic_ref_index >= 1 && ic_ref_index <= length(t)
            curr_time = t[ic_ref_index]
            json_obj["initial_condition_ref_time"] = curr_time
            println("  [Bookkeeping] Stored Initial Condition Ref Time: $curr_time")
        end
    end

    # Run Regression Test
    # Note: Regression test needs t_pivot? It operates on same t, so t_pivot=mean(t) is fine/consistent.
    regression_test_integration(t, tides, ramp_slope, primary_amp, annual_phase, 0.0, 365.242, 0.0, 0.0, impulse_sigma, dc_param, json_obj)

    # Persist bookkeeping changes immediately
    write_params(json_path, json_obj)

    # Use active_coeffs list for reconstruction
    local active_coeffs = Tuple{Int, ComplexF64}[]
    local c_additive = zeros(ComplexF64, 4)
    local M_final
    local best_sec_amp = 0.0
    local best_year_len = 365.242
    local best_primary_amp = primary_amp

    if !skip_optim && staged_path !== nothing
        println("Running Staged Optimization using calibration file: $staged_path")
        backup = archive_params(json_path)
        println("Archived previous params to $backup")
        
        # Read calibration file
        t_cal, I_cal = read_timeseries(staged_path)
        println("Loaded calibration data (n=$(length(t_cal)))")
        
        # Stage 1: Fit Manifold to Calibration File (threshold 0.996)
        kappa_stage1 = collect(-2000:2000)
        
        best_tides, best_ramp, best_primary_amp, best_phi, sec_amp, year_len, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph, best_dc, best_sigma, b_ic_A, b_ic_B, best_r = optimize_manifold_random(t_cal, I_cal, json_obj, ic_A, ic_B, kappa_stage1; manifold_only=true, stop_threshold=0.96, ref_index=ic_ref_index)
        
        println("Stage 1 Complete. Best Calibration Correlation: $best_r")
        
        # Update local vars for Stage 2
        best_sec_amp = sec_amp
        best_year_len = year_len
        
        # Save intermediate params
        json_obj["tides"] = best_tides
        json_obj["ramp_slope"] = best_ramp
        json_obj["primary_impulse_amp"] = best_primary_amp
        json_obj["annual_phase"] = best_phi
        json_obj["secondary_impulse_amp"] = best_sec_amp
        json_obj["year_length"] = best_year_len
        json_obj["drift_linear"] = best_d_lin
        json_obj["drift_quad"] = best_d_quad
        json_obj["comp_annual_amp"] = b_ca_amp
        json_obj["comp_annual_phase"] = b_ca_ph
        json_obj["comp_semi_amp"] = b_csa_amp
        json_obj["comp_semi_phase"] = b_csa_ph
        json_obj["dc_offset"] = best_dc
        json_obj["impulse_sigma"] = best_sigma
        json_obj["initial_condition"] = Dict("A" => b_ic_A, "B" => b_ic_B)
        write_params(json_path, json_obj)
        println("Saved Stage 1 parameters to $json_path")
        
        # Stage 2: Construct Manifold for Target 't' and Run OMP
        println("Stage 2: Applying optimized manifold to target time-series...")
        
        # Construct M_final using clean manifold (inc_add=false)
        inc_add = false 
        M_final = build_manifold(t, best_tides, best_ramp, best_primary_amp, best_phi, sec_amp, year_len, b_ic_A, b_ic_B, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph, best_dc, best_sigma; include_additive=inc_add, ref_index=ic_ref_index)
        M_final .*= 1.0
        
        # Run OMP on target
        kappa = collect(-2000:2000)
        r_final, c_final = run_omp(M_final, I, kappa, t; max_atoms=18)
        
        println("Stage 2 Complete. Final Correlation: $r_final")
        
        # Split coefficients and save
        n_kappa = length(kappa)
        c_manifold = c_final[1:n_kappa]
        c_additive = c_final[n_kappa+1:end]
        
        # Populate active_coeffs for reconstruction
        for (i, val) in enumerate(c_manifold)
            if abs(val) > 1e-9
                push!(active_coeffs, (kappa[i], val))
            end
        end
        
        json_obj["coefficients"] = sparsify_coeffs(kappa, c_manifold)
        json_obj["additive_coefficients"] = Dict(
            "real" => real.(c_additive),
            "imag" => imag.(c_additive)
        )
        write_params(json_path, json_obj)
        println("Saved final coefficients to $json_path")

    elseif !skip_optim && final_mod_path === nothing
        backup = archive_params(json_path)
        println("Archived previous params to $backup")
        
        # Use broader kappa for search
        kappa = collect(-2000:2000)
        
        # Optimize
        best_tides, best_ramp, best_primary_amp, best_phi, sec_amp, year_len, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph, best_dc, best_sigma, b_ic_A, b_ic_B, best_r = optimize_manifold_random(t, I, json_obj, ic_A, ic_B, kappa; manifold_only=manifold_only, ref_index=ic_ref_index, t_pivot=t_pivot_in)
        
        println("Optimization Complete. Best Correlation: $best_r")
        println("Optimized Year Length: $year_len")
        println("Optimized Ramp Slope: $best_ramp")
        println("Optimized Primary Amp: $best_primary_amp")
        println("Optimized Drift Linear: $best_d_lin")
        println("Optimized Drift Quad: $best_d_quad")
        println("Optimized Comp. Annual: Amp=$b_ca_amp, Phase=$b_ca_ph")
        println("Optimized Comp. Semi-Annual: Amp=$b_csa_amp, Phase=$b_csa_ph")
        println("Optimized DC Offset: $best_dc")
        println("Optimized Impulse Sigma: $best_sigma")
        println("Optimized Initial Conditions: A=$b_ic_A, B=$b_ic_B")
        
        best_sec_amp = sec_amp
        best_year_len = year_len
        
        # Reconstruct final manifold using same additive flag as optimization
        inc_add = !manifold_only
        M_final = build_manifold(t, best_tides, best_ramp, best_primary_amp, best_phi, sec_amp, year_len, b_ic_A, b_ic_B, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph, best_dc, best_sigma; include_additive=inc_add, ref_index=ic_ref_index, t_pivot=t_pivot_in)
        M_final .*= 1.0

        # Update JSON (always update manifold params)
        json_obj["tides"] = best_tides
        json_obj["ramp_slope"] = best_ramp
        json_obj["primary_impulse_amp"] = best_primary_amp
        json_obj["annual_phase"] = best_phi
        json_obj["secondary_impulse_amp"] = best_sec_amp
        json_obj["year_length"] = best_year_len
        json_obj["drift_linear"] = best_d_lin
        json_obj["drift_quad"] = best_d_quad
        json_obj["comp_annual_amp"] = b_ca_amp
        json_obj["comp_annual_phase"] = b_ca_ph
        json_obj["comp_semi_amp"] = b_csa_amp
        json_obj["comp_semi_phase"] = b_csa_ph
        json_obj["dc_offset"] = best_dc
        json_obj["impulse_sigma"] = best_sigma
        json_obj["initial_condition"] = Dict("A" => b_ic_A, "B" => b_ic_B)
        
        if manifold_only
             # If manifold only, we stop here (no OMP, no coeff update)
             write_params(json_path, json_obj)
             println("Saved optimized manifold parameters to $json_path")
             
             # Save manifold time series to m.dat (tab-delimited: t  value)
             open("m.dat", "w") do io
                 writedlm(io, [t M_final])
             end
             println("Saved manifold timeseries to m.dat")
             
             println("Manifold Only Complete.")
             return
        end
        
        r_final, c_final = run_omp(M_final, I, kappa, t; max_atoms=18)
        
        # Split coefficients
        n_kappa = length(kappa)
        c_manifold = c_final[1:n_kappa]
        c_additive = c_final[n_kappa+1:end]
        
        # Populate active_coeffs for reconstruction
        for (i, val) in enumerate(c_manifold)
            if abs(val) > 1e-9
                push!(active_coeffs, (kappa[i], val))
            end
        end
        
        json_obj["coefficients"] = sparsify_coeffs(kappa, c_manifold)
        json_obj["additive_coefficients"] = Dict(
            "real" => real.(c_additive),
            "imag" => imag.(c_additive)
        )
        write_params(json_path, json_obj)
        println("Saved optimized parameters to $json_path")
        
    elseif final_mod_path !== nothing
        println("Running Final Modulation Only using manifold from: $final_mod_path")
        
        # Read manifold timeseries
        t_man, M_man = read_timeseries(final_mod_path)
        
        # Alignment
        man_dict = Dict(round(ti, digits=5) => val for (ti, val) in zip(t_man, M_man))
        M_final = zeros(Float64, length(t)) # M_final_aligned replaced by M_final directly
        matched_count = 0
        for (i, ti) in enumerate(t)
            key = round(ti, digits=5)
            if haskey(man_dict, key)
                M_final[i] = man_dict[key]
                matched_count += 1
            else
                M_final[i] = 0.0 # Warning?
            end
        end
        if matched_count < length(t)
            println("Warning: Timeseries mismatched or incomplete overlap ($matched_count / $(length(t))). Interpolated/Zero-filled.")
        end
        
        # Run OMP with expanded range for faster modulation
        println("Using expanded kappa range (-2000:2000) for final modulation.")
        kappa = collect(-2000:2000) 
        r_final, c_final = run_omp(M_final, I, kappa, t; max_atoms=18)
        
        println("Final Modulation Complete. Correlation: $r_final")
        
        # Split coefficients
        n_kappa = length(kappa)
        c_manifold = c_final[1:n_kappa]
        c_additive = c_final[n_kappa+1:end]
        
        # Populate active_coeffs for reconstruction
        for (i, val) in enumerate(c_manifold)
            if abs(val) > 1e-9
                push!(active_coeffs, (kappa[i], val))
            end
        end
        
        json_obj["coefficients"] = sparsify_coeffs(kappa, c_manifold)
        json_obj["additive_coefficients"] = Dict(
            "real" => real.(c_additive),
            "imag" => imag.(c_additive)
        )
        write_params(json_path, json_obj)
        println("Saved final coefficients to $json_path")
        
    else
        println("Skipping optimization, using parameters from $json_path")
        # Read new params if available
        best_sec_amp = Float64(get(json_obj, "secondary_impulse_amp", 0.0))
        best_year_len = Float64(get(json_obj, "year_length", 365.242))
        best_ramp = Float64(get(json_obj, "ramp_slope", 0.001))
        best_primary_amp = Float64(get(json_obj, "primary_impulse_amp", 1.0))
        d_lin = Float64(get(json_obj, "drift_linear", 0.0))
        d_quad = Float64(get(json_obj, "drift_quad", 0.0))
        
        ca_amp = Float64(get(json_obj, "comp_annual_amp", 0.0))
        ca_ph = Float64(get(json_obj, "comp_annual_phase", 0.0))
        csa_amp = Float64(get(json_obj, "comp_semi_amp", 0.0))
        csa_ph = Float64(get(json_obj, "comp_semi_phase", 0.0))
        
        dc_val = Float64(get(json_obj, "dc_offset", 0.0))
        sigma = Float64(get(json_obj, "impulse_sigma", 0.0128))
        
        M_final = build_manifold(t, tides, best_ramp, best_primary_amp, annual_phase, best_sec_amp, best_year_len, ic_A, ic_B, d_lin, d_quad, ca_amp, ca_ph, csa_amp, csa_ph, dc_val, sigma; ref_index=ic_ref_index)
        M_final .*= 1.0
        
        # Load active coeffs
        # In Julia, to assign to outer local active_coeffs, we need no 'local' prefix if already defined
        active_coeffs = load_active_coeffs(json_obj)
        
        if haskey(json_obj, "additive_coefficients")
             add_re = Float64.(json_obj["additive_coefficients"]["real"])
             add_im = Float64.(json_obj["additive_coefficients"]["imag"])
             c_additive = add_re .+ im .* add_im
        else
             c_additive = zeros(ComplexF64, 4)
        end
             
        if isempty(active_coeffs) && !haskey(json_obj, "coefficients")
            println("Warning: No coefficients found in JSON. Running OMP once.")
            _, c_final = run_omp(M_final, I, kappa, t; max_atoms=18)
            n_kappa = length(kappa)
            c_manifold = c_final[1:n_kappa]
            c_additive = c_final[n_kappa+1:end]
            
            # Populate active_coeffs
            for (i, val) in enumerate(c_manifold)
                if abs(val) > 1e-9
                    push!(active_coeffs, (kappa[i], val))
                end
            end
        end
        
        # Check if coefficients are in old dense format (Dict with arrays)
        # If so, convert to new sparse format and save
        coeffs_data = get(json_obj, "coefficients", nothing)
        if coeffs_data isa Dict && haskey(coeffs_data, "kappa")
             println("Converting dense coefficients to sparse format...")
             dense_c = convert(Vector{ComplexF64}, coeffs_data["real"] .+ im .* coeffs_data["imag"])
             json_obj["coefficients"] = sparsify_coeffs(coeffs_data["kappa"], dense_c)
             write_params(json_path, json_obj)
             println("Saved converted parameters to $json_path")
        end
    end

    # ---------- Reconstruct for Plotting ----------
    println("Reconstructing model...")
    
    # Check if 11manifold.dat exists for comparison
    M_empirical = fill(NaN, length(t))
    if isfile("11manifold1.dat")
        t_force, I_force = read_timeseries("11manifold1.dat")
        
        # Create dictionary for lookup with rounded keys to avoid float precision issues
        # Monthly data is ~0.0833, so 4 digits is safe
        force_dict = Dict(round(ti, digits=4) => val for (ti, val) in zip(t_force, I_force))
        
        valid_model = Float64[]
        valid_emp = Float64[]
        
        for (i, ti) in enumerate(t)
            key = round(ti, digits=4)
            if haskey(force_dict, key)
                val = force_dict[key]   # * 10.0
                M_empirical[i] = val
                push!(valid_model, M_final[i])
                push!(valid_emp, val)
            end
        end
        
        if !isempty(valid_emp)
            r_manifold = cor(valid_model, valid_emp)
            println("Correlation between Parametric Manifold and Empirical 11manifold1.dat (n=$(length(valid_emp))): $r_manifold")
        else
             println("Warning: No overlapping time points found with 11manifold1.dat")
        end
    end

    open("model_fit.csv", "w") do io
        println(io, "t,I,I_model,Manifold,Manifold_Empirical")
        
        # Calculate model
        Imodel = zeros(ComplexF64, length(t))
        # Manifold part - Sparse Summation
        for item in active_coeffs
             # active_coeffs is list of (k, val) tuples or objects?
             # load_active_coeffs returns Vector{Tuple{Int, ComplexF64}}
             # sparsify_coeffs returns Dict list.
             # In skip_optim, we use active_coeffs which is Tuple.
             
             # If using sparse JSON, active_coeffs is Tuple list.
             # If using dense JSON, load_active_coeffs converts to Tuple list.
             
             k_val = item[1]
             c_val = item[2]
             
             @. Imodel += c_val * exp(1im * k_val * M_final)
        end
        
        # Additive part
        omega_ann = 2π
        # indices 1..4 corresponding to cos(wt), sin(wt), cos(2wt), sin(2wt)
        if length(c_additive) >= 4
            @. Imodel += c_additive[1] * cos(omega_ann * t)
            @. Imodel += c_additive[2] * sin(omega_ann * t)
            @. Imodel += c_additive[3] * cos(2 * omega_ann * t)
            @. Imodel += c_additive[4] * sin(2 * omega_ann * t)
        end
        
        for i in eachindex(t)
            m_emp_val = length(M_empirical) >= i ? M_empirical[i] : 0.0
            println(io, "$(t[i]),$(I[i]),$(real(Imodel[i])),$(real(M_final[i])),$m_emp_val")
        end
    end
    println("Model fit saved to model_fit.csv")
end


if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: script.jl timeseries.txt params.json [--skip-optim] [--manifold] [--final_mod path_to_manifold]")
    else
        skip_optim = "--skip-optim" in ARGS
        manifold_only = "--manifold" in ARGS
        
        final_mod_path = nothing
        idx_fm = findfirst(x -> x == "--final_mod", ARGS)
        if idx_fm !== nothing && length(ARGS) >= idx_fm + 1
            final_mod_path = ARGS[idx_fm + 1]
        end
        
        staged_path = nothing
        idx_staged = findfirst(x -> x == "--staged", ARGS)
        if idx_staged !== nothing && length(ARGS) >= idx_staged + 1
            staged_path = ARGS[idx_staged + 1]
        end
        
        main(ARGS[1], ARGS[2]; skip_optim=skip_optim, manifold_only=manifold_only, final_mod_path=final_mod_path, staged_path=staged_path)
    end
end
