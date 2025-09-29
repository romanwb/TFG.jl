using TFG, ContinuationSuite, DifferentialEquations, ForwardDiff
use_formal_theme!()
" Duffing Oscillator: ẍ + ξ̃x + x + ϵx³ = 2F₀cos(wt)
    ξ̃ = ϵξ
    F₀ = ϵf₀
"

N, H = 2^6, 2
ξ, ϵ, f₀ = 1, 0.01, 1
ξ̃ = ϵ * ξ
F₀ = ϵ * f₀
E, Eᴴ = fft_matrices(N, H)
x̂₀ = zeros(2H+1)

t = range(0, 2π, length = N + 1)[1:(end - 1)]
f̂₀ = zeros(2H+1)
f̂₀[2] = 2F₀

g(x) = x^3
dg(x) = 3x^2

λ₀ = 0.96
# n = 200 #nº de frecuencias evaluadas en la integracion temporal
# y₀ = [0, 0] # CI del problema temporal

p = DuffingParamsContinuation(N, H, ξ̃, ϵ, E, Eᴴ, x̂₀, f̂₀, g, dg)

cont_pars = ContinuationParameters(λmin = λ₀, λmax = 1.04, Δs = 0.005, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

# la funcion "duffing_continuation" (duffing_oscillator.jl) contiene el sistema de ecs
prob = ContinuationProblem(duffing_continuation, cont_pars, p; autodiff = true)
sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ
amplitudes = [maximum(E * sol.x[:, i]) for i in 1:length(λ_values)]
amplitudes_abs = [maximum(abs.(E * sol.x[:, i])) for i in 1:length(λ_values)]
# @save "scripts/duffing/continuation/continuation.jld2" λ_values amplitudes_abs

# Import time integration (axis X, Y)
# Δω_axis, time_integration = time_integration_values(ξ̃, ϵ, n, λ₀, f̂₀, y₀, g)
using JLD2
@load "data/time_integration/w_var_data.jld2" ω_axis time_integration p
@load "scripts/duffing/continuation/continuation.jld2" λ_values amplitudes_abs

#Import asymptotic solution
#deltaOmega_plus, deltaOmega_minus, Avals = asymptotic_curve(ξ, f̂₀)


# Plot 
using GLMakie, TFG
use_formal_theme!()
begin
    fig = Figure()
    ax = Axis(fig[1, 1],   
    limits=(-0.04, 0.04, 0, 2.2),
    xlabel = L"Δω_{f}", 
    ylabel = L"|x(t)|", 
    # xticklabelfont = "Latin Modern Roman",
    # yticklabelfont = "Latin Modern Roman",
    )
    lin = lines!(ax, λ_values .-1 , amplitudes_abs, color = :blue, label="Continuation Method", linewidth = 2)
    sca = scatter!(ax, ω_axis, time_integration;
         marker = :circle,
         color = :white,  
         strokecolor = :blue,      
         strokewidth = 0.7,           
         markersize = 8,
         label = "Time integration"
    )
    # Legend([lin, sca], ["Continuation", "Time Integration"])
    axislegend(ax; position = :lt, labelsize = 16)  # :rt = right-top
    fig
end

using CairoMakie
save("figures/duffing/continuation/continuation.pdf", fig; backend=CairoMakie)

# mathieu equation resonancia parametrica
