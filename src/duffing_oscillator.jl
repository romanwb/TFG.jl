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

    F₀ = x̂[1] + ϵ * ĝ[1] - f̂[1]
    F = A * x̂[2:end] + ϵ * ĝ[2:end] - f̂[2:end]
    return [F₀; F]
end

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