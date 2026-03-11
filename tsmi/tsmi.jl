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
    kappa = collect(Float64.(j["kappa"]))
    annual_phase = Float64(get(j, "annual_phase", 0.0))
    
    # Read initial condition parameters (default to 0.0)
    ic = get(j, "initial_condition", Dict("A"=>0.0, "B"=>0.0))
    ic_A = Float64(get(ic, "A", 0.0))
    ic_B = Float64(get(ic, "B", 0.0))
    
    return tides, damping, kappa, annual_phase, ic_A, ic_B, j
end

function archive_params(path::String)
    ts = Dates.format(now(), "yyyymmdd_HHMMSS")
    backup = "$(path).bak_$ts"
    cp(path, backup; force=true)
    return backup
end

function write_params(path::String, json_obj)
    open(path, "w") do io
        JSON3.write(io, json_obj; indent=2)
    end
end

# ---------- Physics: manifold construction ----------

function build_tidal_signal(t, tides)
    # t is in years, tides period is in days.
    # Convert t to days for tidal calculation
    # Use consistent floating point type
    T = eltype(t)
    t_days = t .* 365.25
    sig = zeros(T, length(t))
    for td in tides
        P  = Float64(td["period"])
        A  = Float64(td["amplitude"])
        φ  = Float64(td["phase"])
        ω  = 2π / P
        
        # User hint: Strongest tidal period is 27.2122 days.
        # We can emphasize this component if needed, but the amplitude should already reflect importance.
        # However, if the fit is struggling, we might want to boost its effect or ensure it's not drowned out.
        # For now, let's just use the amplitude as given.
        @. sig += A * sin(ω * t_days + φ)
    end
    return sig
end

function build_annual_comb(t, phi)
    # Gaussian comb with period 1 year
    # sigma = 1/12 (approx 1 month)
    σ = 1.0/12.0
    T = promote_type(eltype(t), typeof(phi))
    comb = zeros(T, length(t))
    
    # We assume t is in years.
    # range of years to cover
    # Use simpler iteration to avoid type issues with min/max on Duals if necessary
    # But usually min/max works.
    # To be safe with ForwardDiff, use values or simple loops.
    # But for loop boundaries we need integers.
    # We can use the raw values of t for bounds if t is not being differentiated (t is fixed data).
    # t is fixed data here.
    
    y_min = floor(Int, minimum(t)) - 1
    y_max = ceil(Int, maximum(t)) + 1
    
    for k in y_min:y_max
        # Center of Gaussian for year k
        μ = k + phi
        # vectorized Gaussian addition
        @. comb += exp( - (t - μ)^2 / (2 * σ^2) )
    end
    
    return comb
end

function response_kernel(t, ζ, ω0)
    # damped oscillator impulse response
    T = promote_type(eltype(t), typeof(ζ), typeof(ω0))
    h = zeros(T, length(t))
    for i in eachindex(t)
        τ = t[i]
        if τ < 0
            h[i] = 0.0
        else
            h[i] = exp(-ζ * ω0 * τ) * sin(ω0 * τ)
        end
    end
    return h
end

function build_manifold(t, tides, ζ, ω0, phi, ic_A, ic_B)
    base = build_annual_comb(t, phi) .* build_tidal_signal(t, tides)
    dt = t[2] - t[1]
    # build kernel on same grid (centered at 0..max(t))
    h = response_kernel(t .- t[1], ζ, ω0)
    
    # Forced response (convolution)
    M_forced = conv(base, h) .* dt
    M_forced = M_forced[1:length(t)]
    
    # Free response (initial conditions)
    # Homogeneous solution: exp(-ζωt) * (A*sin(ωt) + B*cos(ωt))
    τ = t .- t[1]
    decay = exp.(-ζ .* ω0 .* τ)
    theta = ω0 .* τ
    
    M_free = decay .* (ic_A .* sin.(theta) .+ ic_B .* cos.(theta))

    return M_forced .+ M_free
end

# ---------- Hoyer metric and gradients ----------

function hoyer(a)
    n = length(a)
    return sqrt(n) * sum(abs.(a)) / norm(a) - 1
end

# ---------- Objective and optimization ----------

function objective(x, I, t, tides, kappa; μ=1e-2, λ_reg=1.0)
    nκ = length(kappa)
    nt = length(I)

    # unpack
    c_re = view(x, 1:nκ)
    c_im = view(x, nκ+1:2nκ)
    
    # parameters are optimized in log-space to ensure positivity
    log_ζ    = x[2nκ+1]
    log_ω0   = x[2nκ+2]
    phi      = x[2nκ+3]
    # ic_A     = x[2nκ+4]
    # ic_B     = x[2nκ+5]
    ic_A = 0.0
    ic_B = 0.0

    # Convert back from log-space
    # Check for Inf before exp to avoid DomainError
    if abs(log_ζ) > 100 || abs(log_ω0) > 100
        return Inf
    end

    T = eltype(x)
    ζ    = exp(log_ζ)
    ω0   = exp(log_ω0)

    c = c_re .+ im .* c_im

    # manifold
    M = build_manifold(t, tides, ζ, ω0, phi, ic_A, ic_B)

    # basis and model
    Φ = [exp.(1im * κ .* M) for κ in kappa]
    Imodel = zeros(Complex{eltype(x)}, nt)
    for (j, φ) in enumerate(Φ)
        Imodel .+= c[j] .* φ
    end
    
    # Regularization to keep damping parameters sane
    # Priors: zeta ~ 0.5, omega0 ~ 2.5
    prior_log_ζ = log(0.5)
    prior_log_ω0 = log(2.5)
    # Also regularize ICs lightly to avoid drift if unconstrained
    # reg_term = λ_reg * ((log_ζ - prior_log_ζ)^2 + (log_ω0 - prior_log_ω0)^2 + 0.1*(ic_A^2 + ic_B^2))
    reg_term = λ_reg * ((log_ζ - prior_log_ζ)^2 + (log_ω0 - prior_log_ω0)^2)

    # objective
    I_real = real.(Imodel)
    # Use Correlation Coefficient as primary metric
    # Manual correlation calculation for ForwardDiff compatibility
    I_mean = mean(I)
    Im_mean = mean(I_real)
    I_centered = I .- I_mean
    Im_centered = I_real .- Im_mean
    
    # Avoid division by zero
    norm_I = norm(I_centered)
    norm_Im = norm(Im_centered)
    
    correlation = if norm_I > 1e-9 && norm_Im > 1e-9
        dot(I_centered, Im_centered) / (norm_I * norm_Im)
    else
        0.0
    end
    
    # Minimize (1 - correlation) strongly
    # Keep MSE for scale, but weight it less? Or rely on correlation.
    # We add MSE to keep coefficients from exploding/vanishing.
    
    J = hoyer(abs.(c)) + μ * sum(abs2.(I_real .- I)) + 10.0 * (1.0 - correlation) + reg_term
    
    return J
end

function optimize_state(I, t, tides, damping, kappa, phi_init, ic_A_init, ic_B_init; μ=1e-2)
    nκ = length(kappa)

    ζ0  = Float64(damping["zeta"])
    ω00 = Float64(damping["omega0"])
    
    # Ensure initial values are positive for log-transform
    # If they are not (e.g. from a bad previous run), reset to defaults
    if ζ0 <= 0
        ζ0 = 0.5
    end
    if ω00 <= 0
        ω00 = 2.5
    end

    # initial manifold and coefficients
    M0 = build_manifold(t, tides, ζ0, ω00, phi_init, ic_A_init, ic_B_init)
    Φ0 = [exp.(1im * κ .* M0) for κ in kappa]
    c0 = [dot(I, φ) for φ in Φ0]

    # x layout: [Re(c), Im(c), log(zeta), log(omega0), phi]
    # We must ensure log arguments are positive. 
    # Use real() to strip potential Complex wrapper if any, though input should be Float64.
    x0 = [real.(c0); imag.(c0); log(ζ0); log(ω00); phi_init]

    # result = optimize(x -> objective(x, I, t, tides, kappa; μ=μ, λ_reg=1.0),
    #                   x0, LBFGS(); autodiff = :forward)
    result = optimize(x -> objective(x, I, t, tides, kappa; μ=μ, λ_reg=1.0),
                      x0, LBFGS(); autodiff = :forward)

    x_opt = Optim.minimizer(result)
    c_re = x_opt[1:nκ]
    c_im = x_opt[nκ+1:2nκ]
    
    log_ζ_opt  = x_opt[2nκ+1]
    log_ω0_opt = x_opt[2nκ+2]
    phi        = x_opt[2nκ+3]
    # ic_A       = x_opt[2nκ+4]
    # ic_B       = x_opt[2nκ+5]
    ic_A = 0.0
    ic_B = 0.0
    
    ζ    = exp(log_ζ_opt)
    ω0   = exp(log_ω0_opt)

    c_opt = c_re .+ im .* c_im
    return c_opt, ζ, ω0, phi, ic_A, ic_B
end

# ---------- Main ----------

function main(ts_path::String, json_path::String; μ=1e-2, skip_optim=false)
    t, I = read_timeseries(ts_path)
    tides, damping, kappa, annual_phase, ic_A, ic_B, json_obj = read_params(json_path)

    if !skip_optim
        backup = archive_params(json_path)
        println("Archived previous params to $backup")

        c_opt, ζ_opt, ω0_opt, phi_opt, ic_A_opt, ic_B_opt = optimize_state(I, t, tides, damping, kappa, annual_phase, ic_A, ic_B; μ=μ)

        # update JSON object
        json_obj["damping"]["zeta"]  = ζ_opt
        json_obj["damping"]["omega0"] = ω0_opt
        json_obj["annual_phase"] = phi_opt
        
        # Store optimized initial conditions
        if !haskey(json_obj, "initial_condition")
            json_obj["initial_condition"] = Dict()
        end
        json_obj["initial_condition"]["A"] = ic_A_opt
        json_obj["initial_condition"]["B"] = ic_B_opt
        
        # optionally store c_opt as magnitude/phase or real/imag
        json_obj["coefficients"] = Dict(
            "kappa" => kappa,
            "real"  => real.(c_opt),
            "imag"  => imag.(c_opt),
        )

        write_params(json_path, json_obj)
        println("Updated params written to $json_path")
    else
        println("Skipping optimization, using parameters from $json_path")
        ζ_opt = Float64(damping["zeta"])
        ω0_opt = Float64(damping["omega0"])
        phi_opt = annual_phase
        ic_A_opt = ic_A
        ic_B_opt = ic_B
        
        # Load c_opt from JSON if available
        if haskey(json_obj, "coefficients")
            c_re = Float64.(json_obj["coefficients"]["real"])
            c_im = Float64.(json_obj["coefficients"]["imag"])
            c_opt = c_re .+ im .* c_im
        else
            # Reconstruct c0 if not in json (unlikely if we are skipping opt)
            M0 = build_manifold(t, tides, ζ_opt, ω0_opt, phi_opt, ic_A_opt, ic_B_opt)
            Φ0 = [exp.(1im * κ .* M0) for κ in kappa]
            c_opt = [dot(I, φ) for φ in Φ0]
        end
    end

    # ---------- Reconstruct model for plotting ----------
    M_opt = build_manifold(t, tides, ζ_opt, ω0_opt, phi_opt, ic_A_opt, ic_B_opt)
    Φ_opt = [exp.(1im * κ .* M_opt) for κ in kappa]
    Imodel = zeros(eltype(c_opt), length(t))
    for (j, φ) in enumerate(Φ_opt)
        Imodel .+= c_opt[j] .* φ
    end
    
    # Save to CSV
    open("model_fit.csv", "w") do io
        println(io, "t,I,I_model")
        for i in eachindex(t)
            println(io, "$(t[i]),$(I[i]),$(real(Imodel[i]))")
        end
    end
    println("Model fit saved to model_fit.csv")
end

# run if called as script
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: script.jl timeseries.txt params.json [--skip-optim]")
    else
        skip_optim = length(ARGS) >= 3 && ARGS[3] == "--skip-optim"
        main(ARGS[1], ARGS[2]; skip_optim=skip_optim)
    end
end
