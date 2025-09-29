using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
use_formal_theme!()

path05 = "scripts/realistic_models/turbine_blade_model/plot_data/EO1/7x7_ff5e-1_mu1e-2.jld2"
path1 = "scripts/realistic_models/turbine_blade_model/plot_data/EO1/7x7_ff1_mu1e-2.jld2"
path2 = "scripts/realistic_models/turbine_blade_model/plot_data/EO1/7x7_ff2_mu1e-2.jld2"
path5 = "scripts/realistic_models/turbine_blade_model/plot_data/EO1/7x7_ff5_mu1e-2.jld2"
Ωff05, Xⱼff05, Aⱼff05, pff05 = load(path05, "Ω", "Xⱼ", "Aⱼ", "p")
Ωff1, Xⱼff1, Aⱼff1, pff1 = load(path1, "Ω", "Xⱼ", "Aⱼ", "p")
Ωff2, Xⱼff2, Aⱼff2, pff2 = load(path2, "Ω", "Xⱼ", "Aⱼ", "p")
Ωff5, Xⱼff5, Aⱼff5, pff5 = load(path5, "Ω", "Xⱼ", "Aⱼ", "p")


    H = pff1.H
    # Xⱼ
    # X₁ = Xⱼ[1:(2H + 1), :]
    # # X₁ = Xⱼ[40*(2H + 1) + 1:41*(2H + 1), :]
    #  X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    # E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    # x₁ = E * X₁
    # x_max = vec(maximum(abs.(x₁), dims = 1))

    #  lines(Ω, x_max)
     A₁ff05 = Aⱼff05[1:(2H + 1), :]
     A₁ff1 = Aⱼff1[1:(2H + 1), :]
     A₁ff2 = Aⱼff2[1:(2H + 1), :]
     A₁ff5 = Aⱼff5[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ff05 = E * A₁ff05
    a₁ff1 = E * A₁ff1
    a₁ff2 = E * A₁ff2
    a₁ff5 = E * A₁ff5

    a_maxff05 = vec(maximum(abs.(a₁ff05), dims = 1))
    a_maxff1 = vec(maximum(abs.(a₁ff1), dims = 1))
    a_maxff2 = vec(maximum(abs.(a₁ff2), dims = 1))
    a_maxff5 = vec(maximum(abs.(a₁ff5), dims = 1))
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, Ωff1, a_maxff1; label=L"1.0\text{f_{ext}}")
    lines!(ax, Ωff05, a_maxff05)
    lines!(ax, Ωff2, a_maxff2)
    lines!(ax, Ωff5, a_maxff5)
    limits!(0.8, 1.05, 0, 110)
    fig
    using CairoMakie
    save("figures/tbm_cyclic/comparison.pdf", fig; backend=CairoMakie)


