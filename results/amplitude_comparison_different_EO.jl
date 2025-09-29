using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
use_formal_theme!()
EO = ("1", "2", "10", "20", "30", "60")
paths = Vector{String}(undef, length(EO))
for i in eachindex(EO)
    paths[i] = "scripts/realistic_models/turbine_blade_model/plot_datav2/EO$(EO[i])/7x7_ff1_mu1e-2.jld2"
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


    H = 5
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
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    # for i in eachindex(paths)
        for i in eachindex(paths)
            lines!(ax1, Ω[i], a_max[i]; label="EO = $(EO[i])")
            lines!(ax2, p[i].ω₀² * Ω[i], p[i].α * a_max[i]; label="EO = $(EO[i])")
        end
    # lines!(ax, Ω2, a_max[i]; label=L"\text{345x345}")
    axislegend()
    # limits!(0.8, 1.05, 0, 110)
    fig
end
    #  using CairoMakie
    #  save("figures/tbm_cyclic/comparison_different_EO.pdf", fig; backend=CairoMakie)

    # =================================================================================================================== #

 
