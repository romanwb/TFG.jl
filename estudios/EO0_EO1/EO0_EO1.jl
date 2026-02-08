using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, LaTeXStrings, Printf
use_formal_theme!()

EOs = [0, 1]


begin
    fig = Figure(resolution = (1000, 800))

    # EO = EOs[n]
    dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
    data_set_selected = 1
    paths = Vector{String}(undef, length(EOs))

    for i in eachindex(EOs)
        EO = EOs[i]
        paths[i] = "scripts/tbm_cyclic_data/all_EO/7x7_EO$(EO).jld2"
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
        ax1 = Axis(fig[1, 1], limits = (nothing, nothing, nothing, nothing), xlabel = "Ω", ylabel = L"||a_{1}||")
        for i in eachindex(paths)
            colors = [:gray, :blue]
            # forcing_float = [0.5, 1, 2, 4]
            # lines!(ax1, Ω[i]*p[i].ω₀², a_max[i]*p[i].α; label = "EO = $(EOs[i])")
            lines!(ax1, Ω[i], a_max[i]; label = "EO = $(EOs[i])", color = colors[i])
        end
        axislegend()
    end
    display(GLMakie.Screen(), fig)
end
limits!(0.96, 1.025, nothing, nothing)
# fig
using CairoMakie
save("estudios/EO0_EO1/EO0_EO1.pdf", fig; backend=CairoMakie)
