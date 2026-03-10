#!/usr/bin/env julia

using DelimitedFiles
using JSON3
using Dates
using LinearAlgebra
using Optim

# ---------- I/O helpers ----------

function read_timeseries(path::String)
    data = readdlm(path)  # whitespace-delimited
    t = data[:, 1]
    I = data[:, 2]
    return t, I
end

function read_params(path::String)
    txt = read(path, String)
    j = JSON3.read(txt)
    tides = j["tides"]          # array of objects with period, amplitude, phase
    damping = j["damping"]      # object with zeta, omega0
    kappa = collect(Float64.(j["kappa"]))
    return tides, damping, kappa, j
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
    sig = zeros(length(t))
    for td in tides
        P  = Float64(td["period"])
        A  = Float64(td["amplitude"])
        φ  = Float64(td["phase"])
        ω  = 2π / P
        sig .+= A .* sin.(ω .* t .+ φ)
    end
    return sig
end

function build_annual_comb(t)
    # simple annual comb: 1 at integer years, 0 otherwise (can refine)
    years = round.(Int, t ./ 365.25)
    comb = zeros(length(t))
    for y in unique(years)
        idx = findall(==(y), years)
        comb[idx] .= 1.0
    end
    return comb
end

function response_kernel(t, ζ, ω0)
    # damped oscillator impulse response
    h = similar(t)
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

function build_manifold(t, tides, ζ, ω0)
    base = build_annual_comb(t) .* build_tidal_signal(t, tides)
    dt = t[2] - t[1]
    # build kernel on same grid (centered at 0..max(t))
    h = response_kernel.(t .- t[1], ζ, ω0)
    M_full = conv(base, h) .* dt
    return M_full[1:length(t)]
end

# ---------- Hoyer metric and gradients ----------

function hoyer(a)
    n = length(a)
    return sqrt(n) * sum(abs.(a)) / norm(a) - 1
end

function hoyer_grad_c(c)
    a = abs.(c)
    n = length(c)
    l1 = sum(a)
    l2 = norm(a)
    dH_da = sqrt(n) * (l2 .* sign.(a) .* l2 .- l1 .* (a ./ l2)) ./ (l2^2)
    return dH_da .* (c ./ (a .+ 1e-12))
end

# ---------- Objective and optimization ----------

function objective_and_grad(x, I, t, tides, kappa; μ=1e-2)
    nκ = length(kappa)
    nt = length(I)

    # unpack
    c_re = view(x, 1:nκ)
    c_im = view(x, nκ+1:2nκ)
    ζ    = x[2nκ+1]
    ω0   = x[2nκ+2]

    c = c_re .+ im .* c_im

    # manifold
    M = build_manifold(t, tides, ζ, ω0)

    # basis and model
    Φ = [exp.(1im * κ .* M) for κ in kappa]
    Imodel = zeros(ComplexF64, nt)
    for (j, φ) in enumerate(Φ)
        Imodel .+= c[j] .* φ
    end

    # objective
    J = hoyer(abs.(c)) + μ * sum(abs2.(Imodel .- I))

    # grad wrt c
    g_c = hoyer_grad_c(c)
    resid = Imodel .- I
    for (j, φ) in enumerate(Φ)
        g_c[j] += 2μ * dot(resid, φ)
    end

    # grad wrt ζ, ω0 via finite differences (simple, can replace with AD)
    ϵ = 1e-5
    # ζ
    M_ζp = build_manifold(t, tides, ζ + ϵ, ω0)
    Φ_ζp = [exp.(1im * κ .* M_ζp) for κ in kappa]
    Imodel_ζp = zeros(ComplexF64, nt)
    for (j, φ) in enumerate(Φ_ζp)
        Imodel_ζp .+= c[j] .* φ
    end
    dJ_dζ = μ * (sum(abs2.(Imodel_ζp .- I)) - sum(abs2.(Imodel .- I))) / ϵ

    # ω0
    M_ωp = build_manifold(t, tides, ζ, ω0 + ϵ)
    Φ_ωp = [exp.(1im * κ .* M_ωp) for κ in kappa]
    Imodel_ωp = zeros(ComplexF64, nt)
    for (j, φ) in enumerate(Φ_ωp)
        Imodel_ωp .+= c[j] .* φ
    end
    dJ_dω0 = μ * (sum(abs2.(Imodel_ωp .- I)) - sum(abs2.(Imodel .- I))) / ϵ

    # pack gradient
    g = [real.(g_c); imag.(g_c); dJ_dζ; dJ_dω0]
    return J, g
end

function optimize_state(I, t, tides, damping, kappa; μ=1e-2)
    nκ = length(kappa)

    ζ0  = Float64(damping["zeta"])
    ω00 = Float64(damping["omega0"])

    # initial manifold and coefficients
    M0 = build_manifold(t, tides, ζ0, ω00)
    Φ0 = [exp.(1im * κ .* M0) for κ in kappa]
    c0 = [dot(I, φ) for φ in Φ0]

    x0 = [real.(c0); imag.(c0); ζ0; ω00]

    result = optimize(x -> objective_and_grad(x, I, t, tides, kappa; μ=μ),
                      x0, LBFGS(); autodiff = :forward)

    x_opt = Optim.minimizer(result)
    c_re = x_opt[1:nκ]
    c_im = x_opt[nκ+1:2nκ]
    ζ    = x_opt[2nκ+1]
    ω0   = x_opt[2nκ+2]

    c_opt = c_re .+ im .* c_im
    return c_opt, ζ, ω0
end

# ---------- Main ----------

function main(ts_path::String, json_path::String; μ=1e-2)
    t, I = read_timeseries(ts_path)
    tides, damping, kappa, json_obj = read_params(json_path)

    backup = archive_params(json_path)
    println("Archived previous params to $backup")

    c_opt, ζ_opt, ω0_opt = optimize_state(I, t, tides, damping, kappa; μ=μ)

    # update JSON object
    json_obj["damping"]["zeta"]  = ζ_opt
    json_obj["damping"]["omega0"] = ω0_opt
    # optionally store c_opt as magnitude/phase or real/imag
    json_obj["coefficients"] = Dict(
        "kappa" => kappa,
        "real"  => real.(c_opt),
        "imag"  => imag.(c_opt),
    )

    write_params(json_path, json_obj)
    println("Updated params written to $json_path")
end

# run if called as script
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: script.jl timeseries.txt params.json")
    else
        main(ARGS[1], ARGS[2])
    end
end
