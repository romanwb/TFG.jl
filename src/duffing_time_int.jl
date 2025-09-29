function time_integration_values(ξ̃, ϵ, n, λ₀, f̂₀, y₀, g, Δω_axis)
    tspan = (0.0, 600.0) # !
    diff, tol = Inf, 1e-2
    time_integration = zeros(n)
    for i ∈ 1:n
        old_max = 0
        p = DuffingParamsTimeIntegration(ξ̃, ϵ, f̂₀, i, n, λ₀, g, Δω_axis)
        tspan = (0.0, 600.0)
        tol = 1e-4
        diff = Inf
        while diff > tol
            step = 25
            prob = ODEProblem(duffing_time_domain, y₀, tspan, p)
            sol = DifferentialEquations.solve(prob, Tsit5(), saveat=0.1) # !!
            new_max = maximum(abs.(sol[1, end-step:end]))
            diff, old_max = abs(new_max - old_max), new_max
            tspan = (tspan[1], tspan[2] + step)
        end
        time_integration[i]=old_max
    end

    # Δω_axis = zeros(n) # !!! 
    # for i in 1:n
    # Δω_axis[i] = 1-(1-λ₀)*cos(0.5*π*(2*(i-1)/(n-1)))
    # end

    return time_integration
end


# ! tspan: intervalo de integracion temporal
# !! saveat: se controla cada cuanto tiempo el integrador devuelve valores
# !!! LinRange es un tipo especial, si no se asigna previamente un vector vacio a la variable, no se genera un vector, se genera un LinRange