" Funciones correspondientes al modelo de fricción"


@inline function normal_force(xₙ, kₙ, xₙ₀)
    u = xₙ - xₙ₀
    if u > 0.0
        return kₙ * u
    else
        return zero(xₙ)
    end
end


@inline function tangencial_force(xₜ, wₜ, kₜ, μ, N)
    if N > 0.0
        T = kₜ * (xₜ - wₜ)
        if abs(T) < μ * N
            #Stick
            return T, wₜ
        else 
            #Slip
            sg = sign(T)
            wₜ = xₜ - sg * μ * N / kₜ
            T = sg * μ * N
            return T, wₜ
        end
    else
        #Lift-off
        wₜ = xₜ
        return zero(xₜ), wₜ
    end 
end


function g_friction(x, kₜ, kₙ, xₙ₀, μ)
    w = zero(x[1]) # Modificado respecto Javier
    t = range(0, 2π, length = length(x) + 1)[1:(end-1)]
    # N = 1.0 .+ 1.25 .* sin.(t) # N(t) ≠ cte
    #N = 1.0 # N(t) = cte
    N = kₙ * xₙ₀
    t_force = similar(x) # Modificado respecto Javier

    for j in 1:2
        for i in eachindex(x)
            t_force[i], w = tangencial_force(x[i], w, kₜ, μ, N)
        end
    end

    return t_force
end


function g_friction_time(x, t)
    N = 1.0 + 1.25 * sin(t)
    μ = 0.25
    kₜ = 1.0
    _, f = tangencial_force(x, 0.0, kₜ, μ, N)
    return f
end