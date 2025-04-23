using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie

" Duffing Oscillator: ẍ + ξ̃x + x + ϵx³ = 2F₀cos(wt)
    ξ̃ = ϵξ
    F₀ = ϵf₀
"


N, H = 2^6, 2
ξ, ϵ = 1, 0.01
ξ̃ = ϵ * ξ
E, Eᴴ = fft_matrices(N, H)
x̂₀ = zeros(2H+1)

t = range(0, 2π, length = N + 1)[1:(end - 1)]
#P = 0.2
#F = @. P * cos(t)
#f̂₀ = Eᴴ * F

        f̂₀ = zeros(2H+1)
        f̂₀[2] = 0.02
#x̂₀ = [0.0, 0.0, 0.0, 0.0, 0.0]
#f̂₀ = [0.0, 0.02, 0.0, 0.0, 0.0]

g(x) = x^3
dg(x) = 3x^2

#g(x) = g_friction(x)
#dg(x) = autodiff_jac(g, x)
#dg(x) = finite_diff_jac(g, x)
λ₀ = 0.96
n = 200 #nº de frecuencias evaluadas en la integracion temporal
y₀ = [0, 0] # CI del problema temporal

# la funcion "duffing_continuation" (duffing_oscillator.jl) contiene el sistema de ecs
p = DuffingParamsContinuation(N, H, ξ̃, ϵ, E, Eᴴ, x̂₀, f̂₀, g, dg)

cont_pars = ContinuationParameters(λmin = λ₀, λmax = 2.0, Δs = 0.01, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)


prob = ContinuationProblem(duffing_continuation, cont_pars, p; autodiff = true)

sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ

#PONER VALOR ABSOLUTO
amplitudes = [maximum(E * sol.x[:, i]) for i in 1:length(λ_values)]
amplitudes_abs = [maximum(abs.(E * sol.x[:, i])) for i in 1:length(λ_values)]


# Import time integration (axis X, Y)
g_time(x, t) = g(x)
Δω_axis, time_integration = time_integration_values(ξ̃, ϵ, n, λ₀, f̂₀, y₀, g)

#Import asymptotic solution
#deltaOmega_plus, deltaOmega_minus, Avals = asymptotic_curve(ξ, f̂₀)



# Plot 
function plot_comparison(y, pos)
    #ω_plus  = 1 .+ 0.01 .* deltaOmega_plus
    #ω_minus = 1 .+ 0.01 .* deltaOmega_minus
    #ω_axis  = 1 .+ 0.01 .* Δω_axis

        ax = Axis(fig[pos, 1], limits=(0.96, 1.03599, 0, 2.25), xlabel=L"\text{Frecuencia}", ylabel= L"x(t)",
        xticklabelfont = "Latin Modern Roman",
        yticklabelfont = "Latin Modern Roman",
         )

        lines!(ax, λ_values, y, color = :dodgerblue, label="Continuation Method")
        #lines!(ax, ω_plus, [2*A for A in Avals], color = :black, label="Asymptotic")
        #lines!(ax, ω_minus, [2*A for A in Avals], color = :black)
        #lines!(ax, ω_axis, time_integration, linestyle = :dash,
        #    marker = :circle,
        #    color = :black,  
        #    label="Time Integration"
        #    strokecolor = :dodgerblue,      
        #    strokewidth = 0.5,           
        #    markersize = 10,
        #    label = "Time integration"
        #)

        axislegend(ax, position = :lt)
    fig
end
fig = Figure()
fig = plot_comparison(amplitudes_abs, 1)
using CairoMakie
TFG.save_figure_pdf("scripts/continuation_timeint_resonance.pdf", figure)
# mathieu equation resonancia parametrica
