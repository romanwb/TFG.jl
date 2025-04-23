using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

function build_A_hbm(M, C, K, H, ω)
    n = size(M, 1)
    Nh = 2H + 1
    T = typeof(ω)
    A = zeros(T, n * Nh, n * Nh)

    A[1:n, 1:n] .= Matrix{T}(I, n, n)  # componente media

    for k in 1:H
        idx_b = n * (2k - 1) + 1
        idx_a = n * (2k) + 1

        Ab = - (k*ω)^2 * M + k*ω * C + K
        Aa = - (k*ω)^2 * M - k*ω * C + K

        A[idx_b:idx_b+n-1, idx_b:idx_b+n-1] .= Ab
        A[idx_a:idx_a+n-1, idx_a:idx_a+n-1] .= Aa
    end

    return A
end


function duffing_continuation_real(x̂, λ, p::DuffingParamsContinuationReal)
    E, Eᴴ, H = p.E, p.Eᴴ, p.H
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    g = p.g

    n = 5
    Nh = 2H + 1
    T = eltype(x̂)

    # Parámetros físicos
    m = T(1.0)
    k = T(10.0)

    # Matrices físicas
    M = Diagonal(fill(m, n)) |> Matrix
    K = zeros(T, n, n)
    for i in 1:n
        K[i, i] += (i in (1, n)) ? k : 2k
        if i < n
            K[i, i+1] -= k
            K[i+1, i] -= k
        end
    end
    C = 0.05 * K

    # No linealidad proyectada (sin funciones auxiliares)
    N = size(E, 1)
    Xreal = zeros(T, N, n)
    for i in 1:n
        Xreal[:, i] .= E * view(x̂, (i-1)*Nh + 1 : i*Nh)
    end

    gX = zeros(T, N, n)
    for j in 1:N
        gX[j, :] .= g(view(Xreal, j, :))
    end

    ĝ = similar(x̂)
    for i in 1:n
        ĝ[(i-1)*Nh + 1 : i*Nh] .= Eᴴ * gX[:, i]
    end

    # Matriz global A
    A = build_A_hbm(M, C, K, H, λ)

    # Residuo completo
    return A * x̂ + ϵ * ĝ - f̂
end





N, H = 2^6, 2
ξ, ϵ = 0.05, 1
ξ̃ = ϵ * ξ
E, Eᴴ = fft_matrices(N, H)
n = 5
Nh = 2H + 1
m = 1.0
    k = 10.0

# Calcula precarga media estática
Fp_static = [1.0, 1.0, 0.0, 0.0, 0.0]  # o la que estés usando
K = zeros(n, n)
for i in 1:n
    K[i, i] += (i in (1, n)) ? k : 2k
    if i < n
        K[i, i+1] -= k
        K[i+1, i] -= k
    end
end

# penalización para evitar singularidad si hace falta
K[1, 1] += 1e-3

x₀ = K \ Fp_static  # solución del sistema estático
           # solución estática
x̂₀ = zeros(n * Nh)
x̂₀[1:n] .= x₀



t = range(0, 2π, length = N + 1)[1:(end - 1)]
P = [1, 1, 0, 0, 0]
F = P * (sin.(t) .+ sin.(2 .* t))'
f̂₀ = vec(Eᴴ * F')  # tamaño 5*(2H+1)

       # f̂₀ = zeros(2H+1)
       # f̂₀[2] = 0.02
#x̂₀ = [0.0, 0.0, 0.0, 0.0, 0.0]
#f̂₀ = [0.0, 0.02, 0.0, 0.0, 0.0]

function g(x)
    T = eltype(x)  # compatible con Dual
    out = zeros(T, 5)
    out[3:5] .= x[3:5].^3
    return out
end

function dg(x)
    T = eltype(x)
    out = zeros(T, 5)
    out[3:5] .= 3 .* x[3:5].^2
    return out
end



#g(x) = g_friction(x)
#dg(x) = autodiff_jac(g, x)
#dg(x) = finite_diff_jac(g, x)
λ₀ = 0.01
n = 200 #nº de frecuencias evaluadas en la integracion temporal
y₀ = [0, 0] # CI del problema temporal

# la funcion "duffing_continuation" (duffing_oscillator.jl) contiene el sistema de ecs
p = DuffingParamsContinuationReal(N, H, ξ̃, ϵ, E, Eᴴ, x̂₀, f̂₀, g, dg)

cont_pars = ContinuationParameters(λmin = λ₀, λmax = 2.0, Δs = 0.01, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)


prob = ContinuationProblem(duffing_continuation_real, cont_pars, p; autodiff = true)

sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ

function get_xt(sol, E, n, H, dof::Int)
    Nh = 2H + 1
    N = size(E, 1)
    [E * sol.x[(dof-1)*Nh+1 : dof*Nh, i] for i in 1:length(sol.λ)]
end

x3_t = get_xt(sol, E, 5, 2, 3)  # n=5, H=2, dof=3

x3_amp = [maximum(abs.(E * sol.x[(3-1)*Nh+1 : 3*Nh, i])) for i in 1:length(sol.λ)]


# Plot 
function plot_comparison(y, pos)
    #ω_plus  = 1 .+ 0.01 .* deltaOmega_plus
    #ω_minus = 1 .+ 0.01 .* deltaOmega_minus
    #ω_axis  = 1 .+ 0.01 .* Δω_axis

        ax = Axis(fig[pos, 1],  xlabel=L"\text{Frecuencia}", ylabel= L"x(t)",
        xticklabelfont = "Latin Modern Roman",
        yticklabelfont = "Latin Modern Roman",
         )

        lines!(ax, λ_values, y, color = :dodgerblue, label="Continuation Method")
        #lines!(ax, ω_plus, [2*A for A in Avals], color = :black, label="Asymptotic")
        #lines!(ax, ω_minus, [2*A for A in Avals], color = :black)
        #lines!(ax, ω_axis, time_integration, linestyle = :dash,
        #    marker = :circle,
        #    color = :black,  
         #   label="Time Integration"
        #    strokecolor = :dodgerblue,      
        #    strokewidth = 0.5,           
        #    markersize = 10,
        #    label = "Time integration"
        

        axislegend(ax, position = :lt)
    fig
end
fig = Figure()
fig = plot_comparison(x3_amp, 1)
using CairoMakie
TFG.save_figure_pdf("scripts/continuation_timeint_resonance.pdf", figure)
# mathieu equation resonancia parametrica
