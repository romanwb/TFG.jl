using TFG, ContinuationSuite, GLMakie, DifferentialEquations

N, H = 2^6, 2
ξ, ϵ = 0.01, 0.01
E, Eᴴ = fft_matrices(N, H)

x̂₀ = [0.0, 0.0, 0.0, 0.0, 0.0]
λ₀ = 0.96
f̂ = [0.0, 0.02, 0.0, 0.0, 0.0]

g(x) = x^3
dg(x) = 3x^2

# la funcion "duffing_continuation" (duffing_oscillator.jl) contiene el sistema de ecs
p = DuffingParamsContinuation(N, H, ξ, ϵ, E, Eᴴ, x̂₀, f̂, g, dg)

cont_pars = ContinuationParameters(λmin = 0.96, λmax = 1.04, Δs = 0.01, maxsteps = 3000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

    prob = ContinuationProblem(
        (x, λ, p) -> duffing_continuation(x, λ, p),  # Aseguramos el tercer argumento
        cont_pars,
        p; 
        autodiff = true
    )
#prob = ContinuationProblem(duffing_continuation, cont_pars; autodiff = true)
#x₀, λ₀ = [0.0, 5.0], 0.96

sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ
amplitudes = [maximum(E * sol.x[:, i]) for i in 1:length(λ_values)]


deltaOmega_plus, deltaOmega_minus, Avals, Δω_axis, old_max_values = DuffingTimeAsymp()
ω_plus  = 1 .+ 0.01 .* deltaOmega_plus
ω_minus = 1 .+ 0.01 .* deltaOmega_minus
ω_axis  = 1 .+ 0.01 .* Δω_axis


fig = Figure()
ax = Axis(fig[1, 1], limits=(0.96, 1.04, 0, 2.25), xlabel="Frecuencia (Ω)", ylabel="Amplitud",
          title="Curva de resonancia")

lines!(ax, λ_values, amplitudes, color = :dodgerblue, label="Continuation Method")
lines!(ax, ω_plus, [2*A for A in Avals], color = :black, label="Asymptotic")
lines!(ax, ω_minus, [2*A for A in Avals], color = :black)
scatter!(ax, ω_axis, old_max_values;
    marker = :circle,
    color = :white,  
    strokecolor = :dodgerblue,      
    strokewidth = 0.5,           
    markersize = 10,
    label = "Time integration"
)

axislegend(ax, position = :lt)
fig
