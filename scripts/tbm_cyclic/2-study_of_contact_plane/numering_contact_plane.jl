# Para obtener la gráfica base...
using GLMakie, TFG
use_formal_theme!()

begin
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1];
        topspinevisible = false,
        rightspinevisible = false,
        leftspinevisible = false,
        bottomspinevisible = false,
    )

    # --- Dibuja los ejes cartesianos ---
    lines!(ax, [0, 7.7], [0, 0], color=:black, linewidth=2)
    lines!(ax, [0, 0], [0, 7.7], color=:black, linewidth=2)

    # --- Añade puntas de flecha como triángulos ---
    tri_x = Point2f[(8, 0), (7.7, 0.1), (7.7, -0.1)]
    tri_y = Point2f[(0, 8), (-0.1, 7.7), (0.1, 7.7)]

    poly!(ax, tri_x, color=:black)
    poly!(ax, tri_y, color=:black)

    hidedecorations!(ax; label=false)
        ax.xlabel = "t₁"
        ax.ylabel = "t₂"
        ax
        ax.xlabelsize = 30
        ax.ylabelsize = 30

        x = 1:7
        y = 1:7
        points = Point2f.((x[i], y[j]) for i in eachindex(x) for j in eachindex(y))
        scatter!(ax, points, strokecolor=:black, strokewidth = 0.5, color=:lightgray, markersize = 15)

    labels = collect(1:49)
    string.(labels)
    text!(ax, points, text = string.(labels), align = (:left, :top), fontsize = 24)

    fig
end

using CairoMakie
save("scripts/tbm_cyclic/1.first_section/3. study_of_contact_plane/numbering_contact_plane.pdf", fig; backend=CairoMakie)
