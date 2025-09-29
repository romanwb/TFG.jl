using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
    use_formal_theme!()
    @load "data/tbm/FEM_tbm_21-7.jld2" contenido # (Info en raw, sin hacer CB, contiene reacciones, etc)
    @load "data/tbm/data_n_modes5/reduced_ns_29_EO1.jld2" Maxᵏ_full Mxaᵏ_full Mxxᵏ_full Kxxᵏ_full fₐᵏ_full fₓᵏ_full ωₐ²ᵏ_full
    sqrt(ωₐ²ᵏ_full[2][1])
     Rₓ = sparsevec(contenido["Rc"])
     for i in 3:3:147
        Rₓ[i] = Rₓ[i] - 100.0
     end
     Rₓ
    # Rₓ = ones(length(fₓᵏ_full[1])) #se obtiene una curva identica a usar la Rₓ original
    # Rₓ = zeros(length(fₓᵏ_full[1]))
    #  Rₓ = vec(contenido["Rc"])
    g = contenido["g"]
    xc = contenido["xc"]

    # N = 29
    # n_modes = 10
    EO = 1
    Na = 5
    Nx = 147
    Nc = 49
    Ns = 29
    # Ns = 8
    Ne = 36771

     kₙ, kₜ₁, kₜ₂, μ = 1e8, 1e6, 2e6, 0.01
    #  kₙ, kₜ₁, kₜ₂, μ = 1e11, 8e10, 2e10, 0.1

    # ξ = 0.035
    ξ = 1e-6
    # γ = 5e20
    γ = 5e12
    H = 10
    N=2^5

# dimensionless strategy
# begin
#     ω₀² = ωₐ²ᵏ_full[EO + 1][1]
#     α = abs(fₐᵏ_full[EO + 1][1] / ω₀²)
#     # β = abs(0.1 * maximum(abs.(Rₓ)) / kₜ₁) # β for coulomb
#     β = α # β for hard spring

#     fₐᵏ_full = fₐᵏ_full / (α * ω₀²)
#     fₓᵏ_full = fₓᵏ_full * β / (α^2 * ω₀²)
#     Rₓ = Rₓ * β / (α^2 * ω₀²)
#     Maxᵏ_full = Maxᵏ_full * (β / α)
#     Mxaᵏ_full = Mxaᵏ_full * (β / α)
#     ωₐ²ᵏ_full = ωₐ²ᵏ_full / ω₀²
#     Mxxᵏ_full = Mxxᵏ_full * (β / α)^2
#     Kxxᵏ_full = Kxxᵏ_full * (β / α)^2 * 1 / ω₀²

#     kₙ, kₜ₁, kₜ₂ = kₙ * (β / α)^2 * 1 / ω₀², kₜ₁ * (β / α)^2 * 1 / ω₀², kₜ₂ * (β / α)^2 * 1 / ω₀²
#     γ = γ * β^4 / (α^2 * ω₀²)
#     ξ = ξ * √(ω₀²)
# end

# ext_force = ComplexHarmonicForcing([(1, ForcingVector(fₐᵏ_full[EO + 1], fₓᵏ_full[EO + 1]))])
ext_force = ComplexHarmonicForcing([(1, ForcingVector(25*(fₐᵏ_full[EO + 1]), 25*(fₓᵏ_full[EO + 1])))])
    # ext_force = HarmonicForcing([(1, ForcingVector(-im / 2 * 0.05fₐ, -im / 2 * 0.05fₓ))])

# cyclic_data = (
#         ωₐ²ᵏ_full = ωₐ²ᵏ_full, Kxxᵏ_full = Kxxᵏ_full, Mxxᵏ_full = Mxxᵏ_full, Maxᵏ_full = Maxᵏ_full,
#         Mxaᵏ_full = Mxaᵏ_full, fₐᵏ_full = fₐᵏ_full, fₓᵏ_full = fₓᵏ_full,
#         Rₓ = Rₓ, Na = Na, Nx = Nx, Nc = Nc, Ns = Ns, Ne = Ne)

cyclic_data = (
        ωₐ²ᵏ_full = ωₐ²ᵏ_full, Kxxᵏ_full = Kxxᵏ_full, Mxxᵏ_full = Mxxᵏ_full, Maxᵏ_full = Maxᵏ_full,
        Mxaᵏ_full = Mxaᵏ_full, fₐᵏ_full = fₐᵏ_full, fₓᵏ_full = fₓᵏ_full,
        Rₓ = Rₓ, Na = Na, Nx = Nx, Nc = Nc, Ns = Ns, Ne = Ne)


cyclic_rom_struct = CyclicROMStructure(cyclic_data, EO, H, ext_force, ξ)


    #  Ω₀, Ωₘₐₓ = sqrt(ωₐ²ᵏ_full[EO + 1][1]) - 100, sqrt(ωₐ²ᵏ_full[EO + 1][1]) + 100
      Ω₀, Ωₘₐₓ = 4950.0, 5100.0
    # Ω₀, Ωₘₐₓ = 0.8, 1.2 #dimensionless
    # Ω₀, Ωₘₐₓ = 2300.0, 2500.0
    # N = 2^9
    N = 2^5
    X_init = zeros((2H + 1) * Nx)



    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 1.0, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-4, ϵₓ = 1e-6), verbose = true, ncols = 5)

    # friction_laws = [Coulomb_2_1D(kₙ = kₙ_vec[i], xₙ₀ = -0.268, kₜ₁ = kₜ₁_vec[i], kₜ₂ = kₜ₂_vec[i], μ = μ)
                    #    for i in 1:Nc]
    friction_laws = [Coulomb_2_1D(kₙ = 1e9, xₙ₀ = 0.0, kₜ₁ = 1e6, kₜ₂ = 2e6, μ = 0.01)
                     for i in 1:Nc]
    #   friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    # hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws)
    hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ
    Aⱼ = compute_cb_amplitudes_cyclic(sol)

     X₁ = Xⱼ[1:(2H + 1), :]
    # X₁ = Xⱼ[40*(2H + 1) + 1:41*(2H + 1), :]
    # X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))

    lines(Ω, x_max)

    A₁ = Aⱼ[1:(2H + 1), :]
    # A₁ = Aⱼ[(2H + 2):(4H + 2), :]
    # A₁ = Aⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    lines(Ω, a_max)
    # lines(Ω, a_max * α)










































function solver_HBM(γ, cyclic_rom_struct, ContPars, HBMCPars, X_init, Ω₀)
    friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ
    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    return a_max, Xⱼ, Ω
end


  begin
    γ_vals = [5e16, 5e17, 5e18, 5e20]
    colors = [:blue, :red, :orange, :green]
    Ω₀, Ωₘₐₓ = 4900.0, 5200.0

    ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.25, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

    fig = Figure()
    ax = Axis(fig[1, 1], yscale = log10, xlabel = L"\Omega", ylabel = L"|a_1(t)|")
    # ax.xticks = ([4900, 5025, 5150], [L"$4900$", L"$5025$", L"$5150$"])
    # ax.yticks = ([10^(-5.5), 10^(-5), 10^(-4.5), 10^(-4)], [L"10^{-5.5}", L"10^{-5}", L"10^{-4.5}", L"10^{-4.5}"])

    for (i, γᵢ) in enumerate(γ_vals)
        ext_force = ComplexHarmonicForcing([(1, ForcingVector(1*(fₐᵏ_full[EO + 1]), 1*(fₓᵏ_full[EO + 1])))])

        cyclic_rom_struct = CyclicROMStructure(cyclic_data, EO, H, ext_force, ξ)

        a_max, Xⱼ, Ω = solver_HBM(γᵢ, cyclic_rom_struct, ContPars, HBMCPars, X_init, Ω₀)

        label_str = L"\gamma = %$(γᵢ)"
        lines!(ax, Ω, a_max; color = colors[i], label = label_str)
    end

    vlines!(ax, [5064.0]; color = :black, linestyle = :dashdot, linewidth = 1)

    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("figures/tbm_cyclic/gamma_var.pdf", fig)
end