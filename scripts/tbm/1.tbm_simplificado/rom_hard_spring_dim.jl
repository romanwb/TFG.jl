using HBMContinuation, LinearAlgebra, Julianim, CairoMakie, DifferentialEquations, TFG

begin
    rom_data = read_rom_data_csv(joinpath(@__DIR__, "../../../data/rom/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    ξ = 1e-6
    γ = 5e16

    ω₀² = ωₐ²[1]
    α = abs(fₐ[1] / ω₀²)
    β = α

    γ = γ * β^4 / (α^2 * ω₀²)
    fₐ = fₐ / (α * ω₀²)
    fₓ = fₓ * β / (α^2 * ω₀²)
    Rₓ = Rₓ * β / (α^2 * ω₀²)
    Max = Max * (β / α)
    ωₐ² = ωₐ² / ω₀²
    ξ = ξ * √(ω₀²)
    Mxx = Mxx * (β / α)^2
    Kxx = Kxx * (β / α)^2 * 1 / ω₀²

    ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)
end

begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 4920.0 / sqrt(ω₀²), 5050.0 / sqrt(ω₀²)
    # Ω₀, Ωₘₐₓ = 4920.0, 5050.0
    ω_axis = Ω₀:1e-3:Ωₘₐₓ
    N, H = 2^5, 2
    X_init = zeros((2H + 1) * Nx)
end
ω_axis = Ω₀:1e-3:Ωₘₐₓ
length(ω_axis)
begin
    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.1, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

    friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = HBMCProblem(rom_str, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ

    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))
    # lines(Ω, a_max * α)
end

X₁ = Xⱼ[(0(2H + 1) + 1):(1(2H + 1)), :]
E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
x₁ = E * X₁
x₁_max = vec(maximum(abs.(x₁), dims = 1))

# Time integration
i = 1


# using Parameters
# @with_kw mutable struct DuffingParamsTimeIntegration

#     ξ::Float64
#     Max::Matrix{Float64}
#     Mxx::Matrix{Float64}
#     Kxx::Matrix{Float64}
#     ωₐ²::Vector{Float64}
#     fₐ::Vector{Float64}
#     fₓ::Vector{Float64}
#     Rₓ::Vector{Float64}
#     ω_axis::Any
#     g::Function
#     i::Int
#     γ::Float64
# end

# function duffing_time_domain(ẏ, y, p::DuffingParamsTimeIntegration, t)
#     ξ, Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ, ω_axis, γ = p.ξ, p.Max, p.Mxx, p.Kxx, p.ωₐ², p.fₐ, p.fₓ, p.Rₓ, p.ω_axis, p.γ
#     i = p.i

#     Id = Matrix{eltype(ωₐ²)}(I, length(ωₐ²), length(ωₐ²))
#     M = [Id Max
#         Max' Mxx]

#     Kaa = Matrix(Diagonal(ωₐ²))
#     Z = zeros(size(Max))
#     K = [Kaa Z
#         Z' Kxx]
    
#     f = [fₐ;fₓ]
#     Ω = ω_axis

#     n = length(y) ÷ 2
#     x  = @view y[1:n]
#     v  = @view y[n+1:end]
#     dx = @view ẏ[1:n]
#     dv = @view ẏ[n+1:end]

#     dx .= v
#     rhs =  f .* sin(Ω[i]*t)
#     rhs .-= ξ .* (K * v)
#     rhs .-= K * x
#     rhs .-= γ .* (x .^ 3)               # g(x) debe devolver Vector del tamaño n

#     dv .= M \ rhs
# end

# function time_integration_values(p::DuffingParamsTimeIntegration, y₀, ω_axis)
#     ξ, Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ, ω_axis, g, γ = p.ξ, p.Max, p.Mxx, p.Kxx, p.ωₐ², p.fₐ, p.fₓ, p.Rₓ, ω_axis, p.g, p.γ
#     tspan = (0.0, 600.0) # !
#     diff, tol = Inf, 1e-2
#     time_integration = zeros(length(ω_axis))
#         i = 1
#         ampl_prev = 0
#         p = DuffingParamsTimeIntegration(ξ, Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ, ω_axis, g, i, γ)
#         t0 = 0.0
#         ΔT = 100*(2π/ω_axis[1])               # o 100·(2π/Ω) si quieres “N períodos”
#         tol = 1e-4
#         diff = Inf
#         T = 2π/ω_axis[1]
#         while diff > tol
#             step = 25
#             prob = ODEProblem((dy, y, p, t) -> duffing_time_domain(dy, y, p, t), y₀, (t0, t0+ΔT), p)
#             sol  = DifferentialEquations.solve(prob, Tsit5(), saveat=T/50)
#             y₀ = sol.u[end]
#             t0 +=  ΔT
#         # mide amplitud en la última ventana (mejor por tiempo, no por índice)
#         t_cut = sol.t[end] - 10.0          # últimos 10 s, por ejemplo
#         sel = findall(t -> t >= t_cut, sol.t)
#         a_new = maximum(abs.(sol[1, sel])) # o norma por DOF si prefieres

#         if abs(a_new - ampl_prev) < tol
#             break
#         end
#         ampl_prev = a_new
#     end


#     # Δω_axis = zeros(n) # !!! 
#     # for i in 1:n
#     # Δω_axis[i] = 1-(1-λ₀)*cos(0.5*π*(2*(i-1)/(n-1)))
#     # end

#     return time_integration
# end


# p = DuffingParamsTimeIntegration(ξ, Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ, ω_axis, g, i, γ)
# n = length(ωₐ²) + size(Mxx,1)
# y₀ = zeros(2n)

# result = time_integration_values(p, y₀, ω_axis)

# --- Definimos la fuerza no lineal g(x) ---
abstract type NonLinearity end

struct Cubic2 <: NonLinearity
    γ::Float64
end

struct Coulomb2 <: NonLinearity
    kₜ::Float64
    kₙ::Float64
    xₙ₀::Float64
    μ::Float64
end

function nonlinear_term3(γ, x)
    T = eltype(x)
    g = similar(x); fill!(g, zero(T))
    g[6] = γ * x[6]^3
    g[7] =  γ * x[7]^3
    g[8] =  γ * x[8]^3
    g[9] =  γ * x[9]^3
    g[10] = γ * x[10]^3
    g[11] = γ * x[11]^3
    g[12] = γ * x[12]^3
    g[13] = γ * x[13]^3
    g[14] = γ * x[14]^3
    g[15] = γ * x[15]^3
    g[16] = γ * x[16]^3
    g[17] = γ * x[17]^3
    return g
end


function g_friction2(x, kₜ, kₙ, xₙ₀, μ)
    w = zero(x) # Modificado respecto Javier
    # t = range(0, 2π, length = length(x) + 1)[1:(end-1)]
    # N = 1.0 .+ 1.25 .* sin.(t) # N(t) ≠ cte
    #N = 1.0 # N(t) = cte
    N = kₙ * (x-xₙ₀)
    t_force = 0.0 # Modificado respecto Javier

    # for j in 1:2
        # for i in eachindex(x)
            t_force, w = tangencial_force(x, w, kₜ, μ, N)
            # t_force, w = tangencial_force(x, w, kₜ, μ, N)
        # end
    # end

    return t_force
end


function nonlinear_term(p::Coulomb2, x)
    g = zeros(length(x))
    g[3] = g_friction2(x[3], 1.0, 1.0, 3.5, 1.0)
    g[4] = g_friction2(x[4], 1.0, 1.0, 3.5, 1.0)
    g[5] = g_friction2(x[5], 1.0, 1.0, 3.5, 1.0)
    return g
end

function nonlinear_term_fin(x)
    g = zeros(length(x))
    g[3] = g_friction2(x[3], 1.0, 1.0, 3.5, 1.0)
    g[4] = g_friction2(x[4], 1.0, 1.0, 3.5, 1.0)
    g[5] = g_friction2(x[5], 1.0, 1.0, 3.5, 1.0)
    return g
end


function system_real_problem_int(ẏ, y, p, t, ω)
    n = size(p.M, 1)
    x = y[1:n]
    v = y[n+1:end]
    f_t = external_force(p.fₐ, p.fₓ, t, ω)
    # g_x = nonlinear_term(p.nl,x)
    g_x = nonlinear_term3(p.γ, x)

    ẏ[1:n] = v
    ẏ[n+1:end] = p.M \ (-(p.K)*x -(p.ξ)*(p.K)*v - g_x + f_t)
    #ẏ[n+1:end] = p.M \ (-p.K*x -p.ξ*p.K*v + f_t)
end

T = 2π/ω_axis[1]
    ΔT = 30*(2π/ω_axis[1]) 

function system_solver(M, K, γ, ξ, fₐ, fₓ, ω_axis, ω; tmax=100000.0, tol=1e-2, idx=5, y₀=nothing)
    n = size(M, 1)
    if y₀ === nothing
        y₀ = zeros(2n)
    end
    tspan = (3000.0, tmax)
    params = SystemParams(M, K, γ, ξ, fₐ, fₓ)
    p = SystemParams(M, K, γ, ξ, fₐ, fₓ)

    old_max, diff = 0.0, Inf
    step = 25.0
    sol = nothing
    t0 = 10000.0
    ΔT = 30*(2π/ω_axis[1]) 
    T = 2π/ω_axis[1]
    while diff > tol
        T = 2π/ω_axis[1]
        prob = ODEProblem((dy, y, p, t) -> system_real_problem_int(dy, y, p, t, ω), y₀, (t0, t0+10*ΔT), params)
        # sol = DifferentialEquations.solve(prob, Tsit5(), saveat=T/50)
        sol = DifferentialEquations.solve(prob, Tsit5(); saveat=0.1)

        x_val = sol[idx, :]
        new_max = maximum(abs.(x_val[end-Int(round(ΔT)):end]))
        #new_max = norm(x_val[end-Int(step/0.1):end], Inf)
        diff = abs(new_max - old_max)
        old_max = new_max
        # tspan = (tspan[1], tspan[2] + step)
        y₀ = sol.u[end]
        t0 += 10*ΔT
    end

    x_vals = sol[idx, :]
    t_vals = sol.t
    n_cycles = (ω .* t_vals) ./ (2π)  # ← vector con número de ciclos

    return x_vals, n_cycles, sol
end

function external_force(fₐ, fₓ, t, ω)
    f_mod = [fₐ;fₓ]
    f = zeros(length(f_mod))
    for i in eachindex(f)
        f[i] = f_mod[i]*sin(ω*t)
    end
    return f
end

function amplitude_from_solution(x_vals::Vector{Float64}; n_puntos_finales::Int = 100)
    N = length(x_vals)
    inicio = max(1, N - n_puntos_finales + 1)
    x_final = x_vals[inicio:end]
    return maximum(abs.(x_final))
end



function resonance_curve(M, K, γ, ξ, ω_axis; dof_idx=5)

    amplitudes = zeros(length(ω_axis))

    y₀ = nothing
    for (i, ω) in enumerate(ω_axis)
        # "_" para ignorar la segunda variable
        x_vals, n_cycles, sol = system_solver(M, K, γ, ξ, fₐ, fₓ, ω_axis, ω, idx=dof_idx, y₀=y₀)
        amplitudes[i] = amplitude_from_solution(x_vals)
        y₀ = sol.u[end]
    end

    return ω_axis, amplitudes
end

γ

    Id = Matrix{eltype(ωₐ²)}(I, length(ωₐ²), length(ωₐ²))
    M = [Id Max
        Max' Mxx]

    Kaa = Matrix(Diagonal(ωₐ²))
    Z = zeros(size(Max))
    K = [Kaa Z
        Z' Kxx]

    struct SystemParams
    M::Matrix{Float64}
    K::Matrix{Float64}
    γ::Float64
    ξ::Float64
    fₐ::Vector{Float64}
    fₓ::Vector{Float64}
end

sol5, solx5 = system_solver(M, K, γ, ξ, fₐ, fₓ, ω_axis, ω_axis[1]; idx=1)
amp = amplitude_from_solution(sol5)
# size(solx5)
ω_axis, amplitudes_int = resonance_curve(M, K, γ, ξ, ω_axis, dof_idx=1)



function make_figure(label::LaTeXString)
    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = label)

    # ax.xticks = ([4950, 5000, 5050], [L"$4950$", L"$5000$", L"$5050$"])
    # ax.yticks = ([2.5e-5, 5e-5, 7.5e-5], [L"$2.5 \times 10^{-5}$", L"$5.0 \times 10^{-5}$", L"$7.5 \times 10^{-5}$"])

    return fig, ax
end

fig, ax = make_figure(L"|a_1(t)|")
lines!(ax, Ω, a_max)
# lines!(ax, ω_axis, amplitudes_int)
scatter!(ax, 0.97, amp)
display(GLMakie.Screen(), fig)
save("figures/rom/hard_spring/amp_a1_adim_wTI.pdf", fig; backend=CairoMakie)

# Ω = 4977 229, Ω = 5000, 321

using DelimitedFiles
M = [I Max;
     Max' Mxx]

K = [diagm(ωₐ²) zeros(5, Nx);
     zeros(Nx, 5) Kxx]

f = [fₐ; fₓ]

xp = hbm_prob.p.preload_state.xₚ
gp = hbm_prob.p.preload_state.gₚ

# writedlm("M.csv", M)
# writedlm("K.csv", K)
# writedlm("f.csv", f)
# writedlm("xp.csv", xp)
# writedlm("gp.csv", gp)

begin
    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1])

    ax.xticks = ([0, π / 2, π, 3π / 2, 2π],
        [L"$0$", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"])

    t = range(0, 2π, length = length(a₁[:, 476]))

    lines!(ax, [0, 2π], [0, 0]; color = :black, linestyle = :dash, linewidth = 1)
    lines!(ax, t, x₁[:, 986]; color = :black, linewidth = 6)

    ax.xlabel = L"$\Omega t$"
    ax.ylabel = L"$x_1(t)$"

    # ylims!(ax, (-0.55, 0.55))

    # save("x_1_5000.png", fig)

    display(GLMakie.Screen(), fig)
end
