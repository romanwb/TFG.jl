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

function system_real_problem_int(ẏ, y, p, t, ω)
    n = size(p.M, 1)
    x = y[1:n]
    v = y[n+1:end]
    f_t = external_force(t, ω)
    g_x = nonlinear_term(x, p.γ)

    ẏ[1:n] = v
    #ẏ[n+1:end] = p.M \ (-p.K*x -p.ξ*p.K*v - g_x + f_t)
    ẏ[n+1:end] = p.M \ (-p.K*x -p.ξ*p.K*v + f_t)
end


function system_solver(M, K, γ, ξ, ω; tmax=1000.0, tol=1e-4, idx=5)
    n = size(M, 1)
    y₀ = zeros(2n)
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
        diff = abs(new_max - old_max)
        old_max = new_max
        tspan = (tspan[1], tspan[2] + step)
    end

    x_vals = sol[idx, :]
    t_vals = sol.t
    n_cycles = (ω .* t_vals) ./ (2π)  # ← vector con número de ciclos

    return x_vals, n_cycles
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

    for (i, ω) in enumerate(ω_axis)
        # "_" para ignorar la segunda variable
        x_vals, _ = system_solver(M, K, γ, ξ, ω, idx=dof_idx)
        amplitudes[i] = amplitude_from_solution(x_vals)
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
n_ω = 200  # Número de puntos
ω_range = (0.0, 2.0)

sol5, solx5 = system_solver(M, K, γ, ξ, 0.45; tmax = 1e3, idx=5)
typeof(solx5)
ω_axis, amplitudes_int = resonance_curve(M, K, γ, ξ, n_ω, ω_range, dof_idx=5)

fig = Figure()
ax = Axis(fig[1, 1],
    #limits = (0, 2, 0, 2),
    xlabel = L"\Omega", ylabel = L"\max(x_3(t))",
    title = L"\text{Respuesta del nodo 3 en función de } \Omega",
    xticklabelfont = "Latin Modern Roman",
    yticklabelfont = "Latin Modern Roman"
)
lines!(ax, solx5, abs.(sol5))
#    marker = :circle,
#    color = :black,
#    strokecolor = :black,
#    strokewidth = 0.5,
#    markersize = 3
#)
fig


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




























using GLMakie

# === Datos simulados ===
N = 10251
t_vec = range(0, 5, length=N)
X = [sol[1,:];
     sol[2,:];
     sol[3,:];
     sol[4,:];
     sol[5,:]]

X = reshape(X, 5, N)  # 5 filas, cada columna es un instante


# === Parámetros de dibujo ===
mass_width = 0.4
mass_height = 0.2
gap = 1.2
base_y = 0.0
x_pos = 0.0
wall_x = 0.6
wall_width = 0.03
triangle_size = 0.08

# === Crear figura y eje ===
fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1]; aspect = DataAspect(), limits = (-1, 1.4, -1, 6))
hidexdecorations!(ax)
hideydecorations!(ax)

# === Dibujar masas (bloques azules) ===
mass_rects = [Rect(x_pos - mass_width/2, base_y + (5 - i) * gap, mass_width, mass_height) for i in 1:5]
mass_obs = [Observable(r) for r in mass_rects]
mass_plots = [poly!(ax, r; color = :skyblue) for r in mass_obs]

# === Dibujar muro entero (una sola vez) ===
wall_height = 5 * gap
wall_rect = Rect(wall_x, base_y - 0.5, wall_width, wall_height + 1.0)
poly!(ax, wall_rect, color = :black)

# === Función para generar muelles con extremos rectos ===
function spring_coords_with_bars(x, y1, y2; turns=6, amp=0.1, bar=0.15)
    n = 2 * turns + 1
    L = abs(y2 - y1)
    body = L - 2bar
    ys = range(y1 + bar, y2 - bar, length=n)
    xs = [x + (-1)^i * amp for i in 1:n]
    xs[1] = x
    xs[end] = x
    xs_full = vcat([x, x], xs, [x, x])
    ys_full = vcat([y1, y1 + bar], ys, [y2 - bar, y2])
    
    return xs_full, ys_full
end

# === Inicializar muelles ===
spring_lines = [Observable(Point2f[]) for _ in 1:5]
spring_plots = [lines!(ax, s; color=:gray, linewidth=2) for s in spring_lines]

# === Inicializar barras y triángulos de fricción ===
bar_lines = [lines!(ax, [Point2f(0, 0), Point2f(0, 0)]; color=:red, linewidth=2) for _ in 1:3]
triangles = [poly!(ax, [Point2f(0, 0), Point2f(0, 0), Point2f(0, 0)]; color=:red) for _ in 1:3]

# === Función de actualización por frame ===
function update_masses!(step)
    # Actualizar masas
    for i in 1:5
        y_offset = X[i, step]
        new_rect = Rect(x_pos - mass_width/2, base_y + (5 - i) * gap + y_offset, mass_width, mass_height)
        mass_obs[i][] = new_rect
    end

    # Actualizar muelles
    for i in 1:4
        y1 = mass_obs[i][].origin[2]
        y2 = mass_obs[i+1][].origin[2] + mass_obs[i+1][].widths[2]
        xs, ys = spring_coords_with_bars(x_pos, y2, y1)
        spring_lines[i][] = Point2f.(xs, ys)
    end
    # Muelle base
    y_top = mass_obs[5][].origin[2]
    y_bottom = y_top - gap
    xs, ys = spring_coords_with_bars(x_pos, y_bottom, y_top)
    spring_lines[5][] = Point2f.(xs, ys)

    # Actualizar barras de contacto y triángulos
    for (i, idx) in enumerate(3:5)
        y_block = mass_obs[idx][].origin[2] + mass_height / 2
        p1 = Point2f(x_pos + mass_width / 2, y_block)
        p2 = Point2f(wall_x, y_block)
        bar_lines[i][1][] = [p1, p2]

        # Triángulo en la punta
        center = p2
        tip = Point2f(center[1] + triangle_size, center[2])
        top = Point2f(center[1], center[2] + triangle_size/2)
        bottom = Point2f(center[1], center[2] - triangle_size/2)
        triangles[i][1][] = [tip, top, bottom]
    end
end


for step in 1:N
        # Actualizar masas
        for i in 1:5
            y_offset = X[i, step]
            new_rect = Rect(x_pos - mass_width/2, base_y + (5 - i) * gap + y_offset, mass_width, mass_height)
            mass_obs[i][] = new_rect
        end
    
        # Actualizar muelles
        for i in 1:4
            y1 = mass_obs[i][].origin[2]
            y2 = mass_obs[i+1][].origin[2] + mass_obs[i+1][].widths[2]
            xs, ys = spring_coords_with_bars(x_pos, y2, y1)
            spring_lines[i][] = Point2f.(xs, ys)
        end
        # Muelle base
        y_top = mass_obs[5][].origin[2]
        y_bottom = y_top - gap
        xs, ys = spring_coords_with_bars(x_pos, y_bottom, y_top)
        spring_lines[5][] = Point2f.(xs, ys)
    
        # Actualizar barras de contacto y triángulos
        for (i, idx) in enumerate(3:5)
            y_block = mass_obs[idx][].origin[2] + mass_height / 2
            p1 = Point2f(x_pos + mass_width / 2, y_block)
            p2 = Point2f(wall_x, y_block)
            bar_lines[i][1][] = [p1, p2]
    
            # Triángulo en la punta
            center = p2
            tip = Point2f(center[1] + triangle_size, center[2])
            top = Point2f(center[1], center[2] + triangle_size/2)
            bottom = Point2f(center[1], center[2] - triangle_size/2)
            triangles[i][1][] = [tip, top, bottom]
        end
end

# === Grabar animación ===
record(fig, "muelles_y_friccion_estetica.mp4", 1:N; framerate = 10) do step
    update_masses!(step)
end