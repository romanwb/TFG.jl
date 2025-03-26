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
            return sg * μ * N, wₜ
        end
    else
        #Lift-off
        wₜ = xₜ
        return zero(xₜ), wₜ
    end 

end


function g_friction(x)
    w = zero(x[1]) # Modificado respecto Javier
    t = range(0, 2π, length = length(x) + 1)[1:(end-1)]
    N = 1.0 .+ 1.25 .* sin.(t) # N(t) ≠ cte
    #N = 1.0 # N(t) = cte
    μ = 0.25
    kₜ = 1.0
    t_force = similar(x) # Modificado respecto Javier

    for j in 1:2
        for i in eachindex(x)
            t_force[i], w = tangencial_force(x[i], w, kₜ, μ, N[i])
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


# RECORDATORIO: comentar sobre el bucle j in 1:2 (para hacer la gráfica habría que guardar informacion desde w=0)
function g_completa(x)
    w = 0.0
    t = range(0, 2π, length = length(x) + 1)[1:(end-1)]
    N = 1.0 .+ 1.25 .* sin.(t) # N(t) ≠ cte
    #N = 1.0 # N(t) = cte
    μ = 0.5
    kₜ = 1.0

    t_full = zeros(2*length(x)) # Recoge historico completo
    t_force = zeros(length(x))

    idx = 1
    for j in 1:2
        for i in eachindex(x)
            t_force[i], w = tangencial_force(x[i], w, kₜ, μ, N[i])
            t_full[idx] = t_force[i]
            idx += 1
        end
    end

    return t_full
end

x_1 = range(0, 1, length=100)
x_2 = range(1, -1, length=100)
x_3 = range(-1, 0, length=100)
x_full = vcat(collect(x_1), collect(x_2), collect(x_3))

T = g_completa(x_full)


fig = Figure()
ax  = Axis(fig[1,1], xlabel="Contact Displacement", ylabel="Friction Force")

lines!(ax, x_full, T[1:300])
lines!(ax, x_full, T[301:600])
fig
