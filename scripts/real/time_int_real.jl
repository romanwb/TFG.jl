using DifferentialEquations, LinearAlgebra, GLMakie

# --- Definición de parámetros ---
struct SystemParams
    M::Matrix{Float64}
    K::Matrix{Float64}
    γ::Float64
    ξ::Float64
end

# --- Definimos la fuerza externa, dependiente de ω ---
function external_force(t, ω)
    return [sin(ω*t) + sin(2*ω*t); sin(ω*t) + sin(2*ω*t); 0.0; 0.0; 0.0]
end

# --- Definimos la fuerza no lineal g(x) ---
function nonlinear_term(x, γ)
    g = zeros(length(x))
    g[3] = γ * x[3]^3
    g[4] = γ * x[4]^3
    g[5] = γ * x[5]^3
    return g
end

# --- Sistema de ecuaciones de primer orden ---
function coupled_system!(ẏ, y, p, t, ω)
    n = size(p.M, 1)
    x = y[1:n]
    v = y[n+1:end]
    f_t = external_force(t, ω)
    g_x = nonlinear_term(x, p.γ)

    ẏ[1:n] = v
    ẏ[n+1:end] = p.M \ (-p.K*x -p.ξ*p.K*v - g_x + f_t)
end

# --- Función principal: integración en rango de frecuencias ---
function time_integration_values_multinode(M, K, γ, n_ω, ω_range; uniform=true)
    Nx = size(M, 1)  # número de nodos
    y₀ = zeros(2Nx)

    tspan = (0.0, 600.0)
    diff, tol = Inf, 1e-2
    amplitudes = zeros(n_ω)
    ω_axis = zeros(n_ω)

    # Crear el eje de frecuencias
    if uniform
        ω_axis .= LinRange(ω_range[1], ω_range[2], n_ω)
    else
        for i in 1:n_ω
            ω_axis[i] = (1-(4e-2)*cos(0.5*π*(2*(i-1)/(n_ω-1))))*(ω_range[2] - ω_range[1]) + ω_range[1]
        end
    end

    params = SystemParams(M, K, γ, ξ)

    # Bucle sobre cada frecuencia ω
    for i in 1:n_ω
        ω = ω_axis[i]
        old_max = 0
        diff = Inf
        tspan = (0.0, 1000.0)

        while diff > tol
            step = 25
            prob = ODEProblem((dy, y, p, t) -> coupled_system!(dy, y, p, t, ω), y₀, tspan, params)
            sol = DifferentialEquations.solve(prob, Tsit5(), saveat=0.1)

            x5 = sol[5, :]
            new_max = maximum(x5[end-step:end])
            diff, old_max = abs(new_max - old_max), new_max
            tspan = (tspan[1], tspan[2] + step)
        end

        amplitudes[i] = old_max
    end

    return ω_axis, amplitudes
end

# --- Definimos matrices M, K y parámetros del problema ---
n = 5
m = 1.0
k = 10.0
γ = 1
ξ = 0.05

M = m * I(n)
K = k * [
    1.0  -1.0   0.0   0.0   0.0
   -1.0   2.0  -1.0   0.0   0.0
    0.0  -1.0   2.0  -1.0   0.0
    0.0   0.0  -1.0   2.0  -1.0
    0.0   0.0   0.0  -1.0   2.0
]

# --- Ejecutamos la integración ---
n_ω = 100  # Número de puntos
ω_range = (0.0, 2.0)

ω_axis, amplitudes_int = time_integration_values_multinode(M, K, γ, n_ω, ω_range)

# --- Plot bonito ---
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
