using TFG, ContinuationSuite, GLMakie, DifferentialEquations

" Duffing Oscillator: ẍ + ξ̃x + x + ϵx³ = 2F₀cos(wt)
    ξ̃ = ϵξ
    F₀ = ϵf₀
"

N, H = 2^6, 2
ξ, ϵ = 1, 0.01
ξ̃ = ϵ * ξ
E, Eᴴ = fft_matrices(N, H)

x̂₀ = [0.0, 0.0, 0.0, 0.0, 0.0]
f̂₀ = [0.0, 0.02, 0.0, 0.0, 0.0]

g(x) = x^3
dg(x) = 3x^2

λ₀ = 0.96
n = 41 #nº de frecuencias evaluadas en la integracion temporal
y₀ = [0, 0] # CI del problema temporal

# la funcion "duffing_continuation" (duffing_oscillator.jl) contiene el sistema de ecs
p = DuffingParamsContinuation(N, H, ξ̃, ϵ, E, Eᴴ, x̂₀, f̂₀, g, dg)

cont_pars = ContinuationParameters(λmin = 0.96, λmax = 1.04, Δs = 0.01, maxsteps = 3000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

# ! Se necesita el uso de una funcion anónima; duffing_continuation tiene parámetros genéricos
#prob = ContinuationProblem((x, λ, p) -> duffing_continuation(x, λ, p), cont_pars, p; autodiff = true)
prob = ContinuationProblem(duffing_continuation, cont_pars, p; autodiff = true)

sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ
amplitudes = [maximum(abs.(E * sol.x[:, i]) for i in 1:length(λ_values))]


# Import time integration (axis X, Y)
Δω_axis, time_integration = time_integration_values(ξ̃, ϵ, n, λ₀, f̂₀, y₀)

#Import asymptotic solution
deltaOmega_plus, deltaOmega_minus, Avals = asymptotic_curve(ξ, f̂₀)

ω_plus  = 1 .+ 0.01 .* deltaOmega_plus
ω_minus = 1 .+ 0.01 .* deltaOmega_minus
ω_axis  = 1 .+ 0.01 .* Δω_axis

# Plot 
function plot_comparison()
    fig = Figure()
        ax = Axis(fig[1, 1], limits=(0.96, 1.04, 0, 2.25), xlabel="Frecuencia (Ω)", ylabel="Amplitud")

        lines!(ax, λ_values, amplitudes, color = :dodgerblue, label="Continuation Method")
        lines!(ax, ω_plus, [2*A for A in Avals], color = :black, label="Asymptotic")
        lines!(ax, ω_minus, [2*A for A in Avals], color = :black)
        scatter!(ax, ω_axis, time_integration;
            marker = :circle,
            color = :white,  
            strokecolor = :dodgerblue,      
            strokewidth = 0.5,           
            markersize = 10,
            label = "Time integration"
        )

        axislegend(ax, position = :lt)
    fig
end

plot_comparison()