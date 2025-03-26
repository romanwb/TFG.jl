function time_integration_values(ξ̃, ϵ, n, λ₀, f̂₀, y₀, g)
    tspan = (0.0, 300.0)
    diff, tol = Inf, 1e-2
    time_integration = zeros(n)
    for i ∈ 1:n
        old_max = 0
        #p = ( ξ̃, ϵ, F₀, 1-(4e-2)*cos(0.5*π*(2*(i-1)/(n-1))) )
        p = DuffingParamsTimeIntegration(ξ̃, ϵ, f̂₀, i, n, λ₀, g)
        tspan = (0.0, 300.0)
        tol = 1e-4
        diff = Inf
        while diff > tol
            step = 25
            prob = ODEProblem(duffing_time_domain_g, y₀, tspan, p)
            sol = DifferentialEquations.solve(prob, Tsit5(), saveat=0.1)
            new_max = maximum(abs.(sol[1, end-step:end]))
            diff, old_max = abs(new_max - old_max), new_max
            tspan = (tspan[1], tspan[2] + step)
        end
        time_integration[i]=old_max
    end

    Δω_axis=zeros(n)
        for i in 1:n
            Δω_axis[i]=((1-(4e-2)*cos(0.5*π*(2*(i-1)/(n-1))))-1)/ϵ
        end

    return Δω_axis, time_integration
end
