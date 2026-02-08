using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, LaTeXStrings, Printf
use_formal_theme!()
# forcing = ("5e-1", "1", "2", "4")
# friction_name = ["1e-2", "25e-3", "5e-2", "1e-1", "15e-2"]
friction_name = ["15e-2", "1e-1", "5e-2", "25e-3", "1e-2"]
EOs=[1, 2, 3, 5, 10, 30]


begin
    fig = Figure(resolution = (1000, 800))

    EO = EOs[6]
    dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
    data_set_selected = 2
    paths = Vector{String}(undef, length(friction_name))

    for i in eachindex(friction_name)
        paths[i] = "scripts/tbm_cyclic/study_of_friction/data/EO$(EO)/$(dataset[data_set_selected])_mu$(friction_name[i])_mod.jld2"
    end

    begin
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

        H = p[1].H
        E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)

        for i in eachindex(paths)
            X[i] = Xⱼ[i][1:(2H + 1), :]
            A[i] = Aⱼ[i][1:(2H + 1), :]
            a[i] = E * A[i]
            x[i] = E * X[i]
            a_max[i] = vec(maximum(abs.(a[i]), dims = 1))
            x_max = vec(maximum(abs.(x[i]), dims = 1))
        end
    end

    begin
        ax1 = Axis(fig[1, 1], limits = (nothing, nothing, nothing, nothing), title = "EO = $(EO)", xlabel = "Ω", ylabel = L"\text{max}|a_{1}|")
        for i in eachindex(paths)
            # mu_float = [0.01, 0.025, 0.05, 0.1, 0.15]
            mu_float = [0.15, 0.1, 0.05, 0.025, 0.01]
            lines!(ax1, Ω[i], a_max[i]; label = "μ = $(mu_float[i])", linewidth=2.0)
        end
        axislegend(position = :lt)
        
    end
    display(GLMakie.Screen(), fig)
end
limits!(0.96, 1.035, nothing, nothing)
using CairoMakie
save("scripts/tbm_cyclic/study_of_friction/figures/comparison_EO$(EO).pdf", fig; backend=CairoMakie)
