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
    
    return tides, damping, kappa, annual_phase, ic_A, ic_B, ramp_slope, j
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

function build_annual_comb(t, phi, sec_amp, year_len, drift_linear, drift_quad)
    # Gaussian comb with period 365.0 days + drifts
    # t is in years (approx 365.25 days)
    # P_comb = 365.0 / 365.25 approx 0.999... in 't' units?
    # No, let's keep it simple: 
    # The 'nominal' year length in the model is year_len.
    # The impulse happens every 365.0 days.
    # So P_comb = 365.0 / year_len (in units of 't')
    
    P_comb = 365.0 / year_len
    
    σ = 1.0/12.0 
    
    comb = zeros(eltype(t), length(t))
    
    # Range of k to cover t
    # Center time for quadratic drift? Let's use mean(t) to keep numbers small
    t_mean = mean(t)
    
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
        
        # Primary Impulse
        @. comb += exp( - (t - t_actual)^2 / (2 * σ^2) )
        
        # Secondary Impulse
        # Assume secondary impulse shifts with the primary
        t_sec = t_actual + 0.5 * P_comb 
        @. comb += sec_amp * exp( - (t - t_sec)^2 / (2 * σ^2) )
    end
    return comb
end

function build_manifold(t, tides, ramp_slope, phi, sec_amp, year_len, ic_A, ic_B, drift_linear, drift_quad, 
                        comp_ann_amp, comp_ann_phase, comp_semi_amp, comp_semi_phase; include_additive=true)
    # Spin-up: Run for 20 years prior to t[1] to eliminate transient
    
    t_start = t[1]
    dt = t[2] - t[1] 
    if dt <= 0; dt = 1.0/365.242; end
    
    t_spinup = range(t_start - 40.0, stop=t_start - dt, step=dt) 
    
    # Base Multiplicative (Impulse * Tides)
    mult_spinup = build_annual_comb(t_spinup, phi, sec_amp, year_len, drift_linear, drift_quad) .* build_tidal_signal(t_spinup, tides, year_len)
    
    # Additive Compensatory Sinusoids
    omega_ann = 2π
    if include_additive
        comp_spinup = @. comp_ann_amp * sin(omega_ann * t_spinup + comp_ann_phase) + 
                         comp_semi_amp * sin(2 * omega_ann * t_spinup + comp_semi_phase)
        base_spinup = mult_spinup .+ comp_spinup
    else
        base_spinup = mult_spinup
    end
    
    val = 0.0 
    for i in eachindex(t_spinup)
        ramp = copysign(ramp_slope, val)
        val = base_spinup[i] + val - ramp
    end
    
    # Main run
    mult = build_annual_comb(t, phi, sec_amp, year_len, drift_linear, drift_quad) .* build_tidal_signal(t, tides, year_len)
    
    if include_additive
        comp = @. comp_ann_amp * sin(omega_ann * t + comp_ann_phase) + 
                  comp_semi_amp * sin(2 * omega_ann * t + comp_semi_phase)
        base = mult .+ comp
    else
        base = mult
    end
    
    M_forced = zeros(eltype(t), length(t))
    current_val = val
    for i in 1:length(t)
        ramp = copysign(ramp_slope, current_val)
        current_val = base[i] + current_val - ramp
        M_forced[i] = current_val
    end
    
    return M_forced
end



# ---------- Sparsity Helpers ----------

# ---------- OMP Core ----------

function run_omp(M, I, kappa, t; max_atoms=12)
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

# ---------- Optimization ----------

function optimize_manifold_random(t, I, json_obj, ic_A, ic_B, kappa; manifold_only=false)
    println("Starting Canonical Manifold Optimization (Random Search - Linear/Quad Drift)...")
    if manifold_only
        println("  (Mode: Manifold-Only Fit)")
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
    current_sec_amp = Float64(get(json_obj, "secondary_impulse_amp", 0.0))
    current_d_lin = Float64(get(json_obj, "drift_linear", 0.0))
    current_d_quad = Float64(get(json_obj, "drift_quad", 0.0))
    current_ca_amp = Float64(get(json_obj, "comp_annual_amp", 0.0))
    current_ca_phase = Float64(get(json_obj, "comp_annual_phase", 0.0))
    current_csa_amp = Float64(get(json_obj, "comp_semi_amp", 0.0))
    current_csa_phase = Float64(get(json_obj, "comp_semi_phase", 0.0))
    
    # Use input initial conditions as starting point
    current_ic_A = ic_A
    current_ic_B = ic_B
    
    best_tides = deepcopy(current_tides)
    best_ramp = current_ramp
    best_phi = current_phi
    best_year_len = current_year_len
    best_sec_amp = current_sec_amp
    best_d_lin = current_d_lin
    best_d_quad = current_d_quad
    best_ca_amp = current_ca_amp
    best_ca_phase = current_ca_phase
    best_csa_amp = current_csa_amp
    best_csa_phase = current_csa_phase
    best_ic_A = current_ic_A
    best_ic_B = current_ic_B
    
    # Objective function wrapper
    sigma_I = std(I_target) # Compute once outside
    function calc_objective(M_trial, I_ref)
        if manifold_only
            m_mean = mean(M_trial)
            
            # Use precomputed means if possible, but I_ref might change in other contexts? 
            # In manifold_only, I_ref is always I_target.
            i_mean = mean(I_ref)
            
            diff_m = M_trial .- m_mean
            diff_i = I_ref .- i_mean
            
            num = sum(diff_m .* diff_i)
            den = sum(diff_m.^2)
            
            slope = den > 1e-12 ? num / den : 0.0
            
            resid = slope .* diff_m .- diff_i
            mse = mean(resid.^2)
            rmse = sqrt(mse)
            
            r = cor(M_trial, I_ref)
            
            nrmse = rmse / (sigma_I + 1e-9)
            
            # Score = Correlation - 0.5 * NRMSE
            score = r - 0.5 * nrmse
            return score
        else
            # Correlation with OMP reconstruction (General Case)
            # If I_ref is forcing, only need 1 atom to correlate well.
            # But generally we want to fit I_ref.
            # Use simpler OMP for speed? 
            r, _ = run_omp(M_trial, I_ref, kappa, t_target; max_atoms=12)
            return r
        end
    end

    # Determine if we include additive terms in manifold construction
    # If manifold_only is true, we EXCLUDE additive terms (ramp only) per requirement
    inc_add = !manifold_only

    # Initial Baseline
    M = build_manifold(t_target, current_tides, current_ramp, current_phi, current_sec_amp, current_year_len, current_ic_A, current_ic_B, current_d_lin, current_d_quad,
                       current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase; include_additive=inc_add)
    r_best = calc_objective(M, I_target)
    
    println("Initial Correlation: $r_best")
    
    # Pre-Optimization: Grid Search for Annual Phase
    # The annual phase (impulse timing) can have local optima due to narrow impulse width.
    # Scan through one full period (~1.0 year) to find the best basin.
    println("Running Grid Search for Annual Phase...")
    best_phase_scan = current_phi
    r_scan_best = r_best
    
    # 24 steps (twice per month)
    scan_steps = range(0.0, stop=1.0, length=25)[1:end-1] 
    for p_val in scan_steps
        M_scan = build_manifold(t_target, current_tides, current_ramp, p_val, current_sec_amp, current_year_len, current_ic_A, current_ic_B, current_d_lin, current_d_quad,
                               current_ca_amp, current_ca_phase, current_csa_amp, current_csa_phase; include_additive=inc_add)
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

    # Loop
    n_iter = 50000 
    T = 0.1 
    alpha = 0.9995
    
    idx_27 = findfirst(x -> abs(Float64(x["period"]) - 27.2122) < 0.1, current_tides)
    
    # Current state
    curr_tides = deepcopy(best_tides)
    curr_ramp = best_ramp
    curr_phi = best_phi # Ensure curr_phi starts at the best found phase
    curr_year_len = best_year_len
    curr_sec_amp = best_sec_amp
    curr_d_lin = best_d_lin
    curr_d_quad = best_d_quad
    curr_ca_amp = best_ca_amp
    curr_ca_phase = best_ca_phase
    curr_csa_amp = best_csa_amp
    curr_csa_phase = best_csa_phase
    curr_ic_A = best_ic_A
    curr_ic_B = best_ic_B
    curr_r = r_best
    
    for i in 1:n_iter
        cand_tides = deepcopy(curr_tides)
        cand_ramp = curr_ramp
        cand_phi = curr_phi
        cand_year_len = curr_year_len
        cand_sec_amp = curr_sec_amp
        cand_d_lin = curr_d_lin
        cand_d_quad = curr_d_quad
        cand_ca_amp = curr_ca_amp
        cand_ca_phase = curr_ca_phase
        cand_csa_amp = curr_csa_amp
        cand_csa_phase = curr_csa_phase
        cand_ic_A = curr_ic_A
        cand_ic_B = curr_ic_B
        
        # Mutate Tides
        for (idx, td) in enumerate(cand_tides)
            if rand() < 0.05 || (idx == idx_27 && rand() < 0.2)
                scale = 1.0 + randn() * 0.05
                td["amplitude"] = Float64(td["amplitude"]) * scale
                shift = randn() * 0.1
                td["phase"] = Float64(td["phase"]) + shift
            end
        end
        
        # Mutate Global Params
        if rand() < 0.2; cand_ramp = abs(cand_ramp * (1.0 + randn() * 0.1)); end
        if rand() < 0.2; cand_phi += randn() * 0.05; end
        if rand() < 0.2; cand_year_len += randn() * 0.002; end
        if rand() < 0.2; cand_sec_amp = clamp(cand_sec_amp + randn() * 0.05, 0.0, 1.0); end
        
        # Mutate Drifts (Very small steps)
        if rand() < 0.3; cand_d_lin += randn() * 1e-6; end
        if rand() < 0.3; cand_d_quad += randn() * 1e-9; end
        
        # Mutate Compensatory Terms (Only if included)
        if inc_add
            if rand() < 0.2
                cand_ca_amp += randn() * 0.05
                cand_ca_phase += randn() * 0.1
            end
            if rand() < 0.2
                cand_csa_amp += randn() * 0.05
                cand_csa_phase += randn() * 0.1
            end
        end

        # Mutate Initial Conditions
        if rand() < 0.2
            cand_ic_A += randn() * 0.05
            cand_ic_B += randn() * 0.05
        end
        
        # Evaluate
        M_cand = build_manifold(t_target, cand_tides, cand_ramp, cand_phi, cand_sec_amp, cand_year_len, cand_ic_A, cand_ic_B, cand_d_lin, cand_d_quad,
                                cand_ca_amp, cand_ca_phase, cand_csa_amp, cand_csa_phase; include_additive=inc_add)
        r_cand = calc_objective(M_cand, I_target)
        
        # Metropolis-Hastings
        delta = r_cand - curr_r
        if delta > 0 || rand() < exp(delta * 10.0 / T) 
            curr_tides = cand_tides
            curr_ramp = cand_ramp
            curr_phi = cand_phi
            curr_year_len = cand_year_len
            curr_sec_amp = cand_sec_amp
            curr_d_lin = cand_d_lin
            curr_d_quad = cand_d_quad
            curr_ca_amp = cand_ca_amp
            curr_ca_phase = cand_ca_phase
            curr_csa_amp = cand_csa_amp
            curr_csa_phase = cand_csa_phase
            curr_ic_A = cand_ic_A
            curr_ic_B = cand_ic_B
            curr_r = r_cand
            
            if r_cand > r_best
                println("Iter $i: New Best r = $r_cand (was $r_best)")
                r_best = r_cand
                best_tides = deepcopy(cand_tides)
                best_ramp = cand_ramp
                best_phi = cand_phi
                best_year_len = cand_year_len
                best_sec_amp = cand_sec_amp
                best_d_lin = cand_d_lin
                best_d_quad = cand_d_quad
                best_ca_amp = cand_ca_amp
                best_ca_phase = cand_ca_phase
                best_csa_amp = cand_csa_amp
                best_csa_phase = cand_csa_phase
                best_ic_A = cand_ic_A
                best_ic_B = cand_ic_B
            end
        end
        
        T *= alpha
        if i % 500 == 0
            println("  Iter $i... Temp: $(round(T, digits=4)), Curr: $(round(curr_r, digits=4)), Best: $(round(r_best, digits=4))")
        end
    end
    
    return best_tides, best_ramp, best_phi, best_sec_amp, best_year_len, best_d_lin, best_d_quad, best_ca_amp, best_ca_phase, best_csa_amp, best_csa_phase, best_ic_A, best_ic_B, r_best
end

# ---------- Main ----------

function main(ts_path::String, json_path::String; skip_optim=false, manifold_only=false, final_mod_path::Union{String, Nothing}=nothing)
    t, I = read_timeseries(ts_path)
    tides, damping, kappa, annual_phase, ic_A, ic_B, ramp_slope, json_obj = read_params(json_path)

    local c_final
    local M_final
    local best_sec_amp = 0.0
    local best_year_len = 365.242

    if !skip_optim && final_mod_path === nothing
        backup = archive_params(json_path)
        println("Archived previous params to $backup")
        
        # Use broader kappa for search
        kappa = collect(-2000:2000)
        
        # Optimize
        best_tides, best_ramp, best_phi, sec_amp, year_len, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph, b_ic_A, b_ic_B, best_r = optimize_manifold_random(t, I, json_obj, ic_A, ic_B, kappa; manifold_only=manifold_only)
        
        println("Optimization Complete. Best Correlation: $best_r")
        println("Optimized Year Length: $year_len")
        println("Optimized Ramp Slope: $best_ramp")
        println("Optimized Drift Linear: $best_d_lin")
        println("Optimized Drift Quad: $best_d_quad")
        println("Optimized Comp. Annual: Amp=$b_ca_amp, Phase=$b_ca_ph")
        println("Optimized Comp. Semi-Annual: Amp=$b_csa_amp, Phase=$b_csa_ph")
        println("Optimized Initial Conditions: A=$b_ic_A, B=$b_ic_B")
        
        best_sec_amp = sec_amp
        best_year_len = year_len
        
        # Reconstruct final manifold using same additive flag as optimization
        inc_add = !manifold_only
        M_final = build_manifold(t, best_tides, best_ramp, best_phi, sec_amp, year_len, b_ic_A, b_ic_B, best_d_lin, best_d_quad, b_ca_amp, b_ca_ph, b_csa_amp, b_csa_ph; include_additive=inc_add)
        M_final .*= 1.0

        # Update JSON (always update manifold params)
        json_obj["tides"] = best_tides
        json_obj["ramp_slope"] = best_ramp
        json_obj["annual_phase"] = best_phi
        json_obj["secondary_impulse_amp"] = best_sec_amp
        json_obj["year_length"] = best_year_len
        json_obj["drift_linear"] = best_d_lin
        json_obj["drift_quad"] = best_d_quad
        json_obj["comp_annual_amp"] = b_ca_amp
        json_obj["comp_annual_phase"] = b_ca_ph
        json_obj["comp_semi_amp"] = b_csa_amp
        json_obj["comp_semi_phase"] = b_csa_ph
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
        
        r_final, c_final = run_omp(M_final, I, kappa, t; max_atoms=12)
        
        # Split coefficients
        n_kappa = length(kappa)
        c_manifold = c_final[1:n_kappa]
        c_additive = c_final[n_kappa+1:end]
        
        
        json_obj["coefficients"] = Dict(
            "kappa" => kappa,
            "real"  => real.(c_manifold),
            "imag"  => imag.(c_manifold),
        )
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
        r_final, c_final = run_omp(M_final, I, kappa, t; max_atoms=12)
        
        println("Final Modulation Complete. Correlation: $r_final")
        
        # Split coefficients
        n_kappa = length(kappa)
        c_manifold = c_final[1:n_kappa]
        c_additive = c_final[n_kappa+1:end]
        
        json_obj["coefficients"] = Dict(
            "kappa" => kappa,
            "real"  => real.(c_manifold),
            "imag"  => imag.(c_manifold),
        )
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
        d_lin = Float64(get(json_obj, "drift_linear", 0.0))
        d_quad = Float64(get(json_obj, "drift_quad", 0.0))
        
        ca_amp = Float64(get(json_obj, "comp_annual_amp", 0.0))
        ca_ph = Float64(get(json_obj, "comp_annual_phase", 0.0))
        csa_amp = Float64(get(json_obj, "comp_semi_amp", 0.0))
        csa_ph = Float64(get(json_obj, "comp_semi_phase", 0.0))
        
        M_final = build_manifold(t, tides, best_ramp, annual_phase, best_sec_amp, best_year_len, ic_A, ic_B, d_lin, d_quad, ca_amp, ca_ph, csa_amp, csa_ph)
        M_final .*= 1.0
        
        local c_manifold
        local c_additive
        
        if haskey(json_obj, "coefficients")
            c_re = Float64.(json_obj["coefficients"]["real"])
            c_im = Float64.(json_obj["coefficients"]["imag"])
            c_manifold = c_re .+ im .* c_im
            
            if haskey(json_obj, "additive_coefficients")
                add_re = Float64.(json_obj["additive_coefficients"]["real"])
                add_im = Float64.(json_obj["additive_coefficients"]["imag"])
                c_additive = add_re .+ im .* add_im
            else
                c_additive = zeros(ComplexF64, 4)
            end
        else
            println("Warning: No coefficients found in JSON. Running OMP once.")
            _, c_final = run_omp(M_final, I, kappa, t; max_atoms=12)
            n_kappa = length(kappa)
            c_manifold = c_final[1:n_kappa]
            c_additive = c_final[n_kappa+1:end]
        end
    end

    # ---------- Reconstruct for Plotting ----------
    println("Reconstructing model...")
    
    # Check if 11manifold.dat exists for comparison
    M_empirical = fill(NaN, length(t))
    if isfile("11manifold.dat")
        t_force, I_force = read_timeseries("11manifold.dat")
        
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
            println("Correlation between Parametric Manifold and Empirical 11manifold.dat (n=$(length(valid_emp))): $r_manifold")
        else
             println("Warning: No overlapping time points found with 11manifold.dat")
        end
    end

    open("model_fit.csv", "w") do io
        println(io, "t,I,I_model,Manifold,Manifold_Empirical")
        
        # Calculate model
        Imodel = zeros(ComplexF64, length(t))
        # Manifold part
        for (j, k_val) in enumerate(kappa)
             @. Imodel += c_manifold[j] * exp(1im * k_val * M_final)
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
        
        main(ARGS[1], ARGS[2]; skip_optim=skip_optim, manifold_only=manifold_only, final_mod_path=final_mod_path)
    end
end
