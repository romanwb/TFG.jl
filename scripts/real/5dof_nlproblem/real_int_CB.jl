"!
ESTE SCRIPT RESUELVE EL PROBLEMA UNA VEZ SE HA HECHO LA TRANSFORMACION DE CB
"

using DifferentialEquations, LinearAlgebra, GLMakie, Statistics

# --- Definición de parámetros ---
struct SystemParams
    M::Matrix{Float64}
    K::Matrix{Float64}
    γ::Float64
    ξ::Float64
end

# --- Definimos la fuerza externa, dependiente de ω ---
function external_force(t, ω)
    return [(-1.3763819204711734)*(sin(ω*t) + sin(2*ω*t)); (0.3249196962329064)*(sin(ω*t) + sin(2*ω*t)); 3.0*(sin(ω*t) + sin(2*ω*t)); 0.0; 0.0]

end

# --- Definimos la fuerza no lineal g(x) ---
function nonlinear_term(x, γ)
    g = zeros(length(x))
    g[3] = γ * x[3]^3
    g[4] = γ * x[4]^3
    g[5] = γ * x[5]^3
    return g
end

function system_real_problem_int(ẏ, y, p, t, ω)
    n = size(p.M, 1)
    x = y[1:n]
    v = y[n+1:end]
    f_t = external_force(t, ω)
    g_x = nonlinear_term(x, p.γ)

    ẏ[1:n] = v
    ẏ[n+1:end] = p.M \ (-(p.K)*x -(p.ξ)*(p.K)*v - g_x + f_t)
    #ẏ[n+1:end] = p.M \ (-p.K*x -p.ξ*p.K*v + f_t)
end


function system_solver(M, K, γ, ξ, ω; tmax=10000.0, tol=1e-4, idx=5, y₀=nothing)
    n = size(M, 1)
    if y₀ === nothing
        y₀ = zeros(2n)
    end
    tspan = (0.0, tmax)
    params = SystemParams(M, K, γ, ξ)

    old_max, diff = 0.0, Inf
    step = 25.0
    sol = nothing
    while diff > tol
        prob = ODEProblem((dy, y, p, t) -> system_real_problem_int(dy, y, p, t, ω), y₀, tspan, params)
        sol = DifferentialEquations.solve(prob, Tsit5(), saveat=1)

        x_val = sol[idx, :]
        new_max = maximum(abs.(x_val[end-Int(step/0.1):end]))
        #new_max = norm(x_val[end-Int(step/0.1):end], Inf)
        diff = abs(new_max - old_max)
        old_max = new_max
        tspan = (tspan[1], tspan[2] + step)
    end

    x_vals = sol[idx, :]
    t_vals = sol.t
    n_cycles = (ω .* t_vals) ./ (2π)  # ← vector con número de ciclos

    return x_vals, n_cycles, sol
end


function amplitude_from_solution(x_vals::Vector{Float64}; n_puntos_finales::Int = 100)
    N = length(x_vals)
    inicio = max(1, N - n_puntos_finales + 1)
    x_final = x_vals[inicio:end]
    return maximum(abs.(x_final))
end



function resonance_curve(M, K, γ, ξ, n_ω, ω_range; dof_idx=5)

    ω_axis = LinRange(ω_range[1], ω_range[2], n_ω)
    amplitudes = zeros(n_ω)

    y₀ = nothing
    for (i, ω) in enumerate(ω_axis)
        # "_" para ignorar la segunda variable
        x_vals, n_cycles, sol = system_solver(M, K, γ, ξ, ω, idx=dof_idx, y₀=y₀)
        amplitudes[i] = amplitude_from_solution(x_vals)
        y₀ = sol.u[end]
    end

    return ω_axis, amplitudes
end

# --- Definimos matrices M, K y parámetros del problema ---
n = 5
m = 1.0
k = 10.0
γ = 1
ξ = 0.05

omega_a2 = [3.8196601125010514, 26.18033988749895]
#Kaa = Diagonal(omega_a2.^2)
Kaa = Diagonal(omega_a2)

Maa = Matrix(I, 2, 2)

Mxx = [3.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0]

Max = [-1.3763819204711734 0.0 0.0;
        0.3249196962329064 0.0 0.0]

Kxx = [10.0 -10.0 0.0;
       -10.0 20.0 -10.0;
         0.0 -10.0 20.0]

Psi = [1.0 0.0 0.0 1.0 0.0 0.0]

Phi = [-0.850651 -0.525731;
       -0.525731 0.850651]

fa = [-1.3763819204711734, 0.3249196962329064]
fx = [3.0, 0.0, 0.0]
Rx = [-1.2, -0.1, -0.1]
xe0 = zeros(2)

M = [Maa Max;
     Max' Mxx]
K = [Kaa zeros(2,3);
     zeros(3,2) Kxx]
f = zeros(5)
f[1:2] = fa
f[3:5] = fx
f
# --- Ejecutamos la integración ---
n_ω = 100  # Número de puntos
ω_range = (0.0, 2.5)

sol5, solx5 = system_solver(M, K, γ, ξ, 0.63; idx=5)
size(solx5)

ω_axis, amplitudes_int = resonance_curve(M, K, γ, ξ, n_ω, ω_range, dof_idx=5)

begin
fig = Figure()
ax = Axis(fig[1, 1],
    #limits = (0, 2, 0, 2),
    xlabel = L"\Omega", ylabel = L"\max(x_3(t))",
    title = L"\text{Respuesta del nodo 3 en función de } \Omega",
    xticklabelfont = "Latin Modern Roman",
    yticklabelfont = "Latin Modern Roman"
)
lines!(ax, solx5, sol5)
lines!(ax, t_hbm, x_hbm)
#    marker = :circle,
#    color = :black,
#    strokecolor = :black,
#    strokewidth = 0.5,
#    markersize = 3
#)
fig
end

maximum(abs.(sol5))
maximum(abs.(x_hbm))


fig = Figure()
ax = Axis(fig[1, 1],
    limits = (0, 2, 0, 2),
    xlabel = L"\Omega", ylabel = L"\max(x_3(t))",
    title = L"\text{Respuesta del nodo 3 en función de } \Omega",
    xticklabelfont = "Latin Modern Roman",
    yticklabelfont = "Latin Modern Roman"
)
lines!(ax, ω_axis, amplitudes_int)
#    marker = :circle,
#    color = :black,
#    strokecolor = :black,
#    strokewidth = 0.5,
#    markersize = 3
#)
fig
