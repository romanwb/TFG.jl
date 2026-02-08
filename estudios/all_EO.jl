using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, LaTeXStrings, Printf, Colors
use_formal_theme!()
    EOs=Vector{Int}(undef, 74)
    for i in 1:74
        EOs[i]= i-1
    end

    # EOs = [0, 1, 2, 3, 10, 15, 17, 20, 21, 22, 23, 30, 40, 50, 73]

colors = ["#4876b4", "#b44848", "#a3b448","#9648b4", "#4fb448", "#48b4a1",  "#746e4c"]
palette = parse.(Colorant, colors)

ranges = [(1,4), (5,13), (14,30), (31, 74)]
color_of(i) = palette[findfirst(r -> r[1] <= i <= r[2], ranges)]
color_of(2)
begin
    fig = Figure(resolution = (1000, 800))

    # EO = EOs[n]
    dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
    data_set_selected = 1
    paths = Vector{String}(undef, length(EOs))

    for i in eachindex(EOs)
        EO = EOs[i]
        larger_EOs = [0, 1, 2, 3]
        # paths[i] = EO in larger_EOs ? "scripts/tbm_cyclic_data/all_EO/7x7_EO$(EO)_largerOm.jld2" : "scripts/tbm_cyclic_data/all_EO/7x7_EO$(EO).jld2"
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
        ax1 = Axis(fig[1:10, 1:10], limits = (nothing, nothing, nothing, nothing), xlabel = "Ω", ylabel = L"|a_{1}|")
        for i in eachindex(paths)
            # forcing_float = [0.5, 1, 2, 4]
            lines!(ax1, (sqrt.(Ω[i]*p[i].ω₀²))/(2π), a_max[i]*p[i].α, color=color_of(i))
            # text!(ax1, (Ω[i]*p[i].ω₀²)[Int(round(length(Ω[i])/2))], (a_max[i]*p[i].α)[Int(round(length(Ω[i])/2))], text = "EO = $(EOs[i])", align = (:left, :bottom), fontsize = 24)
            # text!(ax1, (Ω[i]*p[i].ω₀², a_max[i]*p[i].α), text = "EO =", align = (:left, :bottom), fontsize = 24)
            # lines!(ax1, Ω[i], a_max[i]; label = "EO = $(EOs[i])")
        end
        # Elementos de leyenda: uno por rango (color representativo)
        labels = ["EO ∈ [$(first(r)),$(last(r))]" for r in ranges]
        elements = [LineElement(color=c, linewidth=3) for c in palette]
        elements_adj = elements[1:length(labels)]
        Legend(fig[1,7], elements_adj, labels; tellwidth=false)

        # axislegend()
    end
    display(GLMakie.Screen(), fig)
end
limits!(0.96, 1.025, nothing, nothing)
# fig
using CairoMakie
save("scripts/tbm_cyclic/study_of_ext_force/figures/comparison_EO$(EO).pdf", fig; backend=CairoMakie)
