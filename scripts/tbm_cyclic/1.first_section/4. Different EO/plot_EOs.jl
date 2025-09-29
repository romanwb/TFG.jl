using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
use_formal_theme!()
EO = ("1", "5", "10", "20", "30", "60")
paths = Vector{String}(undef, length(EO))
dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
key = 9
# for i in eachindex(EO)
#     paths[i] = "scripts/tbm_cyclic/1.first_section/4. Different EO/data/EO$(EO[i])/$(dataset[key]).jld2"
# end

for i in eachindex(EO)
    paths[i] = "scripts/tbm_cyclic/1.first_section/5_different_dataset/data/EO$(EO[i])/$(dataset[key]).jld2"
end

Ω  = Vector{Vector{Float64}}(undef, length(paths))
Xⱼ = Vector{Matrix{Float64}}(undef, length(paths))
X = Vector{Matrix{Float64}}(undef, length(paths))
Aⱼ = Vector{Matrix{Float64}}(undef, length(paths))
A = Vector{Matrix{Float64}}(undef, length(paths))
a = Vector{Matrix{Float64}}(undef, length(paths))
x = Vector{Matrix{Float64}}(undef, length(paths))
a_max = Vector{Vector{Float64}}(undef, length(paths))
p = Vector{Any}(undef, length(paths))

for i in eachindex(paths)
    Ω[i], Xⱼ[i], Aⱼ[i], p[i] = load(paths[i], "Ω", "Xⱼ", "Aⱼ", "p")
end


    H = p[3].H
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    # Xⱼ
    # X₁ = Xⱼ[1:(2H + 1), :]
    # # X₁ = Xⱼ[40*(2H + 1) + 1:41*(2H + 1), :]
    #  X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    # E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    # x₁ = E * X₁
    # x_max = vec(maximum(abs.(x₁), dims = 1))

for i in eachindex(paths)
    X[i] = Xⱼ[i][1:(2H + 1), :]
    A[i] = Aⱼ[i][1:(2H + 1), :]
    a[i] = E * A[i]
    x[i] = E * X[i]
    a_max[i] = vec(maximum(abs.(a[i]), dims = 1))
end

begin
    fig = Figure(resolution = (1200, 800))
    ax1 = Axis(fig[1, 1], xlabel = "Ω [Hz]", ylabel = L"\text{max}|a_{1}|")
    # ax2 = Axis(fig[1, 2])
    # for i in eachindex(paths)
        for i in eachindex(paths)
            # lines!(ax1, Ω[i], a_max[i]; label="EO = $(EO[i])")
            lines!(ax1, (sqrt.(p[i].ω₀² * Ω[i])/2π), p[i].α * a_max[i]; label="EO = $(EO[i])")
            # lines!(ax1, p[i].ω₀² * Ω[i],  a_max[i]; label="EO = $(EO[i])")
        end
    # lines!(ax, Ω2, a_max[i]; label=L"\text{345x345}")
    axislegend(ax1)
    # limits!(0.8, 1.05, 0, 110)
    display(GLMakie.Screen(), fig)
end
    #  using CairoMakie
    #  save("figures/tbm_cyclic/comparison_different_EO.pdf", fig; backend=CairoMakie)

    # =================================================================================================================== #

 limits!(0.90, 1.05, nothing, nothing)
 using CairoMakie
 save("scripts/tbm_cyclic/1.first_section/4. Different EO/comparison_different_EO_$(dataset[key]).pdf", fig; backend=CairoMakie)