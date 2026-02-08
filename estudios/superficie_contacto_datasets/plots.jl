# Para obtener la gráfica base...
using GLMakie, TFG, Colors
use_formal_theme!()
greem = "#48b4a7"
palette = parse.(Colorant, greem)
function print_dataset_contacts(dataset)
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1])

    hidedecorations!()
    #gray points
        x = 1:7
        y = 1:7
        points = Point2f.((x[i], y[j]) for i in eachindex(x) for j in eachindex(y))
        scatter!(ax, points, strokecolor=:black, strokewidth = 0.5, color=:lightgray, markersize = 15)

    #black points
        xs ,ys = dataset_contacts[dataset]
        points = Point2f.((xs[i], ys[j]) for i in eachindex(xs) for j in eachindex(ys))
        scatter!(ax, points, strokecolor=:orange, strokewidth = 2.5, color=:black, markersize = 15)
        fig
return fig
end


# using CairoMakie
# save("scripts/tbm_cyclic/1.first_section/3. study_of_contact_plane/numbering_contact_plane.pdf", fig; backend=CairoMakie)

dataset_contacts = Dict{String, Tuple}()
dataset_contacts["7x7"] = (1:7, 1:7)
dataset_contacts["345x345"] = (3:5, 3:5)
dataset_contacts["23456x23456"] = (2:6, 2:6)
dataset_contacts["34x34"] = (3:4, 3:4)
dataset_contacts["1-7x1-7"] = (1:6:7, 1:6:7)
dataset_contacts["first_last_cols"] = (1:6:7, 1:7)
dataset_contacts["first_last_rows"] = (1:7, 1:6:7)
dataset_contacts["first2colums"] = (1:2, 1:7)
dataset_contacts["upper2rows"] = (1:7, 1:2)


print_dataset_contacts("upper2rows")





dataset_contacts