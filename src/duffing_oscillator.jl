"""
    Duffing oscillator en forma matricial y matriz jacobiana  
"""

# Duffing Oscillator: ẍ + ξẋ + x + ϵg(x) = f(t)
function duffing(x̂, p::DuffingParams)
    E, Eᴴ, H = p.E, p.Eᴴ, p.H
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    g, Ω = p.g, p.Ω

    ĝ = Eᴴ * g.(E * x̂)

    A = system_matrix(H, ξ, Ω)

    #F₀ = x̂[1] + ϵ * ĝ[1] - f̂[1]
    F₀ = x̂[1] + ĝ[1] - f̂[1]
    #F = A * x̂[2:end] + ϵ * ĝ[2:end] - f̂[2:end]
    F = A * x̂[2:end] + ĝ[2:end] - f̂[2:end]
    return [F₀; F]
end

function duffing_continuation(x̂, λ, p::DuffingParamsContinuation)
    E, Eᴴ, H = p.E, p.Eᴴ, p.H
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    g = p.g

    #ĝ = Eᴴ * g.(E * x̂)
    ĝ = Eᴴ * g(E * x̂)
    A = system_matrix(H, ξ, λ)

    #F₀ = x̂[1] + ϵ * ĝ[1] - f̂[1]
    F₀ = x̂[1] + ĝ[1] - f̂[1]
    #F = A * x̂[2:end] + ϵ * ĝ[2:end] - f̂[2:end]
    F = A * x̂[2:end] + ĝ[2:end] - f̂[2:end]
    return [F₀; F]
end

function duffing_time_domain(ẏ, y, p::DuffingParamsTimeIntegration, t)
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    float_f₀ = only(filter(!iszero, f̂))
    i, n, λ₀ = p.i, p.n, p.λ₀
    Ω = 1-(1-λ₀)*cos(0.5*π*(2*(i-1)/(n-1)))

    ẏ[1] = y[2] 
    ẏ[2] = - ξ*y[2] - y[1] - ϵ*y[1]^3 + float_f₀*cos(Ω*t)
end

function duffing_time_domain_g(ẏ, y, p::DuffingParamsTimeIntegration, t)
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    float_f₀ = only(filter(!iszero, f̂))
    i, n, λ₀ = p.i, p.n, p.λ₀
    g = p.g
    Ω = 1-(1-λ₀)*cos(0.5*π*(2*(i-1)/(n-1)))

    ẏ[1] = y[2] 
    ẏ[2] = - ξ*y[2] - y[1] - g(y[1], t) + float_f₀*cos(Ω*t)
end

"
f(du, u, p, t, a, b, c)
f_ode = (du, u, p, t) -> f(du, u, p, t, 0.1, 0.2, 0.3)
"

function jacobian_matrix(x, p::DuffingParams)
    E, Eᴴ, H = p.E, p.Eᴴ, p.H
    ξ, ϵ, f̂ = p.ξ, p.ϵ, p.f̂
    dg, Ω = p.dg, p.Ω

   # Parte lineal: 
    A = system_matrix(H, ξ, Ω)
    J = zeros(2H+1, 2H+1)
    J[1,1] = 1.0
    J[2:end, 2:end] = A

    # Parte no lineal:
    v = E * x
    dv = dg.(v)
    D = Diagonal(dv)
    J_nonlin =  ϵ*(Eᴴ * D * E)

    J .+= J_nonlin
    return J
end