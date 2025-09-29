using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, LaTeXStrings, Printf
use_formal_theme!()
# forcing = ("5e-1", "1", "2", "4")
forcing = ("5e-1", "1", "2", "4")
EOs=[1, 10, 20, 30, 60]

fig = Figure(resolution = (1000, 800))

begin
    EO = EOs[1]
    dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
    data_set_selected = 1
    paths = Vector{String}(undef, length(forcing))

    for i in eachindex(forcing)
        paths[i] = "scripts/tbm_cyclic/study_of_ext_force/data/EO$(EO)/$(dataset[data_set_selected])_ff$(forcing[i]).jld2"
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

        H = 3
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
            forcing_float = [0.5, 1, 2, 4]
            lines!(ax1, Ω[i]*p[1].ω₀², a_max[i]; label = "fₑₓₜ = $(forcing_float[i])")
        end
        axislegend()
    end
    display(GLMakie.Screen(), fig)
end

fig
#using CairoMakie
#save("figures/tbm_cyclic/comparison_different_dataset.pdf", fig; backend=CairoMakie)
