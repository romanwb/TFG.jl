using TFG, ContinuationSuite, GLMakie

function f(x, λ, p)
    x₁, x₂ = x
    f₁ = x₂ - (x₁ + 5)        # y = x + 5 (recta)
    f₂ = x₁^2 + x₂^2 - λ      # x^2 + y^2 = λ (circunferencia)
    return [f₁, f₂]
end

cont_pars = ContinuationParameters(λmin = 5.0, λmax = 20.0, Δs = 0.05, maxsteps = 3000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

prob = ContinuationProblem(f, cont_pars; autodiff = true)
x₀, λ₀ = [0.0, 5.0], 12.6

sol = continuation(prob, x₀, λ₀)
x, λ = sol.x, sol.λ
