using GLMakie
using GeometryBasics: Point2f, widths
const Point2f0 = Point2f
const Vec2f0   = Vec2f

"""
Create the main panel in the multi-panel figure to study the displacements of the contact nodes.

# Arguments
- `Ω`: 1D array of excitation frequencies.
- `a_max`: array storing the maximum displacement amplitude at each value of `Ω`;
  only `a_max[1]` is plotted as the main curve.
- `idxsΩ`: integer indices selecting the frequencies of interest in `Ω` where
  vertical lines and markers are drawn.
- `limits`: 4-element tuple or vector `(xmin, xmax, ymin, ymax)` specifying the
  axis limits.

# Keyword Arguments
- `resolution`: tuple `(width, height)` specifying the figure resolution in pixels.
- `xlabel`: label for the horizontal axis.
- `ylabel`: label for the vertical axis (LaTeXStrings allowed).
- `Ω_labels`: optional vector of custom tick labels associated with `idxsΩ`;
  if `nothing`, the default numeric tick labels are used.
"""
function base_figure(Ω, a_max, idxsΩ, limits;
        resolution = (1400, 1000),
        figure_padding = 120,
        xlabel      = "Ω",
        ylabel      = L"|a_{1}|",
        Ω_labels    = nothing,
    )

    fig = Figure(resolution = resolution, figure_padding=figure_padding, backgroundcolor = :white)

    ax = Axis(fig[2:3,2:3];
        limits = limits,
        xlabel = xlabel,
        ylabel = ylabel,
    )

    # curva principal
    lines!(ax, Ω, a_max[1])

    # puntos de estudio + líneas verticales
    xs = Ω[idxsΩ]
    ys = a_max[1][idxsΩ]
    ymin = limits[3]

    for (x, y) in zip(xs, ys)
        lines!(ax, [x, x], [y, ymin];
            linestyle = :dash, linewidth = 1, color = :gray)
    end

    scatter!(ax, xs, ys; color = :white)

    # ticks
    if Ω_labels !== nothing
        ax.xticks = (xs, Ω_labels)
    end
    ax.xticksmirrored   = false
    ax.xminorgridvisible = false
    ax.xticklabelrotation = π/3
    ax.xticklabelsize     = 20

    return fig, ax
end


"""
Need of info here!!
"""
# ! Solo contempla el caso de 7x7, pendiente de generalizarse
function full_mesh(ax)
    x = 1:7
    y = 1:7
    points = Point2f.((x[i], y[j]) for i in eachindex(x) for j in eachindex(y))
    scatter!(ax, points, strokecolor=:black, strokewidth = 0.5, color=:lightgray, markersize = 15)
        Δ = 0.5
    xlims!(ax, minimum(x) - Δ, maximum(x) + Δ)
    ylims!(ax, minimum(y) - Δ, maximum(y) + Δ)

    hidedecorations!(ax)
    return nothing
end


"""
Need of info here!1
"""
const dataset_distros = Dict{String, Tuple{Vector{Int}, Vector{Int}, Int}}(
    "7x7" => (collect(1:7), collect(1:7), 7),
    "34x34" => ([3, 4], [3, 4], 2),
    "345x345" => ([3, 4, 5], [3, 4, 5], 3),
    "23456x23456" => ([2, 3, 4, 5, 6], [2, 3, 4, 5, 6], 5),
    "1-7x1-7" => ([1,7], [1,7], 2),
    "first_last_cols" => ([1,7], [1,2,3,4,5,6,7], 7),
    "first_last_rows" => ([1,2,3,4,5,6,7], [1,7], 2),
    "first2colums" => ([1,2,3,4,5,6,7], [1,2], 7),
    "upper2rows" => ([1,2], [1,2,3,4,5,6,7], 2),
)

function get_positions(name::AbstractString)
    v1, v2, dim = dataset_distros[name]
    return v1, v2, dim
end

function coords_from_index(n::Int, dim_plot::Int)
    q, r = divrem(n - 1, dim_plot)  # q = cociente, r = resto
    i = q + 1
    j = r + 1
    return i, j
end

# dim refiere a la dimension del cuadrado, que si son Nc=49 es dim=7
function plot_slip_stick(ax, isstick, n, name::AbstractString)
    x, y, dim = get_positions(name)
    isstick_n = isstick[:,n]
    # dim = 2 # "numero de filas"
    dim_plot = 7
    for (i,iₓ) in zip(eachindex(x),x)
        for (j,jₓ) in zip(eachindex(y),y)
            a, b = coords_from_index(Int(dim_plot*(iₓ-1)+jₓ), 7)
            if isstick_n[dim*(i-1)+j] == Int(0)
                scatter!(ax, (a,b), strokecolor=:black, strokewidth = 0.5, color=:black, markersize = 15)
            elseif isstick_n[dim*(i-1)+j] == Int(1)
                scatter!(ax, (a,b), strokecolor=:black, strokewidth = 0.5, color=:LIGHTBLUE, markersize = 15)
            elseif isstick_n[dim*(i-1)+j] == Int(2)
                scatter!(ax, (a,b), strokecolor=:black, strokewidth = 0.5, color=:red, markersize = 15)
            end
        end
    end
end

"""
Need of info here!1
"""
function global_coords_from_local(Nc::Int, name::AbstractString)
    x, y, dim = get_positions(name)        # x,y: índices globales; dim: tamaño del set local

    @assert 1 <= Nc <= dim^2 "Nc = $Nc fuera de rango para dim = $dim"

    # índice local dentro del set (1..dim en cada dirección)
    i_local, j_local = coords_from_index(Nc, dim)

    # índices globales: filas/columnas en la malla 7x7
    i_global = x[i_local]
    j_global = y[j_local]

    return i_global, j_global
end

function plot_selected_contact(ax, Ωs_to_eval, Nc::Int, name::AbstractString)
    x, y = global_coords_from_local(Nc::Int, name::AbstractString)
    for i in eachindex(Ωs_to_eval)
    scatter!(ax, (x,y), strokecolor=:green, strokewidth = 3, color=:transparent, markersize = 22)
    end
end


"De (x,y) en datos de un Axis a un punto en la figura (píxeles)"
function point_in_fig_from_axis(ax::Axis, x::Real, y::Real)
    v = ax.scene.viewport[]
    p_scene = Makie.project(ax.scene, Point2f0(x, y))  # coords en la escena del axis
    # viewport.origin es un Vec, así que lo sumamos y convertimos a Point2f0
    p = p_scene + Vec2f0(v.origin)
    return Point2f0(p[1], p[2])
end

"Centro de un Axis en coords de figura (píxeles)"
function center_of_axis_in_fig(ax::Axis)
    v = ax.scene.viewport[]
    o = v.origin
    w = widths(v)
    return Point2f0(o[1] + w[1]/2,
                    o[2] + w[2]/2)
end




    # row = [1, 1, 1, 4, 4, 4]
    # col = [1, 2, 4, 1, 2, 4]

    # sizing = 0.7
    # fwidth = sizing*ones(6)
    # fheight = sizing*ones(6)
    # fhalign = [0.5, 1.5, 1, 0.5, 1.5, 1]
    # fvalign = [-0.25, 0.5, -0.25, 1.25, 0.5, 1.25]