# Mapa de modos para diferentes set de datos
using TFG, GLMakie, CairoMakie, JLD2
# n_dataset = 9
# ω_set = Vector{Vector{Vector{Float64}}}(undef, n_dataset)
file = "data/tbm/ns146_reduced/reduced_ns146_7x7.jld2"

ω² = load(file, "ωₐ²ᵏ_full")

# dataset_name = ["7x7", "34x34", "345x345", "23456x23456", "first2colums", "upper2rows", "1-7x1-7", "first_last_cols", "first_last_rows"]
x = 1:Int(length(ω²)/2)
y = [sqrt(ω²[i][j]) / 2π for i in x, j in 1:3]

use_formal_theme!() # !
mode = 1
begin
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1]; xlabel="Engine Order", ylabel="Frequency [Hz]")
        for j in 1:3 
            # lines!(ax, x, y[i,:,j], label = "Mode $j, dataset $i")
            scatter!(ax, x, y[:,j], label = "Mode $(j)")
        end
    axislegend(ax; position = :lt)
    display(GLMakie.Screen(), fig)
end

save("scripts/tbm/3.turbine_blade_model/1.first_section/3modes.pdf", fig; backend=CairoMakie)

# ! solo se representa la mitad de los EO ya que se presenta una simetria de espejo. Mirar apuntes
