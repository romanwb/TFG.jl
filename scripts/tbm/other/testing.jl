using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays

# @load "scripts/realistic_models/turbine_blade_model/plot_datav2/EO$(EO)/$(dataset[key])_ff1_mu1e-2.jld2" Ω Xⱼ Aⱼ p

@load "scripts/realistic_models/turbine_blade_model/study_of_ext_force/data/EO1/7x7_ff5e-1.jld2" Ω Xⱼ Aⱼ p
H=5
p
     X₁ = Xⱼ[1:(2H + 1), :]
    # X₁ = Xⱼ[40*(2H + 1) + 1:41*(2H + 1), :]
    # X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^16, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))

    

    A₁ = Aⱼ[1:(2H + 1), :]
    # A₁ = Aⱼ[(2H + 2):(4H + 2), :]
    # A₁ = Aⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^16, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    begin
    f = Figure()
    ax1 = Axis(f[1, 1], limits = (0.91, 0.94, 0, 670), xlabel = "x label", ylabel = "y label",
    title = "Amplitud primer modo, adimensional")
#     ax2 = Axis(f[1, 2], xlabel = "x label", ylabel = "y label",
#     title = "Amplitud primer modo")
#     ax3 = Axis(f[2, 1], xlabel = "x label", ylabel = "y label",
#     title = "Amplitud primer contacto, adimensional")
#     ax4 = Axis(f[2, 2], xlabel = "x label", ylabel = "y label",
#     title = "Amplitud primer contacto")

    lines!(ax1, Ω, a_max)
#     lines!(ax2, Ω * p.ω₀², a_max * p.α)
#     lines!(ax3, Ω, x_max)
#     lines!(ax4, Ω * p.ω₀², x_max * p.α)

    f
    end









    