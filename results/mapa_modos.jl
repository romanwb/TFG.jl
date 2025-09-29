# Mapa de modos para diferentes set de datos
using TFG, GLMakie, CairoMakie, JLD2
n_dataset = 9
ω_set = Vector{Vector{Vector{Float64}}}(undef, n_dataset)
for (i, file) in enumerate([
    "data/tbm/ns146_reduced/reduced_ns146_7x7.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_34x34.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_345x345.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_23456x23456.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_first2colums.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_upper2rows.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_1-7x1-7.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_first_last_cols.jld2"
    "data/tbm/ns146_reduced/reduced_ns146_first_last_rows.jld2"
])
    ω_set[i] = load(file, "ωₐ²ᵏ_full")
end

dataset_name = ["7x7", "34x34", "345x345", "23456x23456", "first2colums", "upper2rows", "1-7x1-7", "first_last_cols", "first_last_rows"]
x = 1:Int(length(ω_set[1])/2)
y = [sqrt(ω_set[i][j][k]) / 2π for i in 1:n_dataset, j in x, k in 1:3]

use_formal_theme!() # !
mode = 1
begin
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1]; xlabel="Engine Order", ylabel="Frequency [Hz]", title="Mode $(mode)")
    for i in 1:n_dataset
        # for j in 1:3 
            # lines!(ax, x, y[i,:,j], label = "Mode $j, dataset $i")
            lines!(ax, x, y[i,:,mode], label = "Dataset $(dataset_name[i])")
        # end
    end
    axislegend(ax; position = :lt)
    fig
end

save("figures/results/mode_shape_mode$(mode).pdf", fig; backend=CairoMakie)

# ! solo se representa la mitad de los EO ya que se presenta una simetria de espejo. Mirar apuntes
