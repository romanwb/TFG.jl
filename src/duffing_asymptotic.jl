function asymptotic_curve(ξ, f̂₀)

    f₀ = 100*(only(filter(!iszero, f̂₀))/2)

    Avals = range(1e-3, 2.2, length=3000)
    deltaOmega_plus  = Float64[]
    deltaOmega_minus = Float64[]

    for A in Avals
        radicando = (f₀/(2*A))^2 - (ξ/2)^2
        if radicando >= 0
            raiz = sqrt(radicando)
            push!(deltaOmega_plus,  (3/2)*A^2 + raiz)
            push!(deltaOmega_minus, (3/2)*A^2 - raiz)
        else
            push!(deltaOmega_plus,  NaN)
            push!(deltaOmega_minus, NaN)
        end
    end

    return deltaOmega_plus, deltaOmega_minus, Avals
end
