using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, ProgressMeter
    use_formal_theme!()

    dataset=("7x7", "345x345", "23456x23456", "34x34", "1-7x1-7", "first_last_cols", "first_last_rows", "first2colums", "upper2rows")
    EOs=[1, 10, 20, 30, 60]
    
# @showprogress 1 "Procesando..."  for m in eachindex(EOs)

        key = 2
        friction_values = [0.01, 0.025, 0.05, 0.1, 0.15]
# for n in eachindex(friction_values)
    @load "data/tbm/ns146_reduced/reduced_ns146_$(dataset[key]).jld2" Maxᵏ_full Mxaᵏ_full Mxxᵏ_full Kxxᵏ_full fₐᵏ_full fₓᵏ_full ωₐ²ᵏ_full Φᵏ Ψᵏ Rₓ
    @load "data/tbm/ns146_raw/ns146_$(dataset[key]).jld2" contenido

        # EO = EOs[m]
        EO = 60
        Ns = 146
        Na = length(Maxᵏ_full[1][:,1])
        Nx = Int(length(Mxxᵏ_full[1][:,1]))
        Nc = Int(Nx/3)
        Ne = size(contenido["Mee_new"], 1)
        # Ne = 39609 # 7x7
        # Ne = 39729 #345x345
        # Ne = 39681 #2345x2345

        # Construccion de las fuerzas
        fe = contenido["F_vector"]
        fe = fe[1:Ne,1]
        fwave = ComplexF64.(fe)  
        Φm, Ψm = Φᵏ[EO+1], Ψᵏ[EO+1]
        fₐ = Φm' * fwave
        fₓ = Ψm' * fwave

        fₓᵏ_full = [zeros(ComplexF64, Nx) for _ in 1:Ns]  # Vector{Vector{ComplexF64}}
        fₓᵏ_full[EO+1] = fₓ
        fₐᵏ_full = [zeros(ComplexF64, Na) for _ in 1:Ns]  # Vector{Vector{ComplexF64}}
        fₐᵏ_full[EO+1] = fₐ


        #   Rₓ = ones(Nx) * -100
        # Construccion de las reacciones
        extra_react = 200.0
        added_force = extra_react*49/Nc

        for i in 3:3:Nx
            # Rₓ[i-2] = Rₓ[i-2] - added_force/extra_react
            # Rₓ[i-1] = Rₓ[i-1] - added_force/extra_react
             Rₓ[i] = Rₓ[i] - added_force
            # Rₓ[i] = Rₓ[i] + 10.0
        end

for i in 1:Int(length(Rₓ))
    println("coefficient: ", i, "value: ", Rₓ[i])
end


    begin
        kr = 49/Nc # parámetro de ajuste
        kₙ, kₜ₁, kₜ₂, μ = 1e9*kr, 1e6*kr, 1e6*kr, friction_values[n]
        # kₙ, kₜ₁, kₜ₂, μ = 1e9*kr, 1e6*kr, 2e6*kr, 0.01
        # kₙ, kₜ₁, kₜ₂, μ = 1e11, 8e10, 2e10, 0.1
        # kₙ, kₜ₁, kₜ₂, μ = 148871.069, 148871.069, 148871.069, 0.1

        # ξ = 0.035
        ξ = 1e-6
        # γ = 5e16
         γ = 5e12

        H = 5 # numero de armonicos
    end

    # adimensionalizacion
    begin
        # ω₀² = 2.0456228983408358e7
        ω₀² = ωₐ²ᵏ_full[EO + 1][1]
        α = abs(fₐᵏ_full[EO + 1][1] / ω₀²)
        β = abs(0.1 * maximum(abs.(Rₓ)) / kₜ₁)
        # β = α

        fₐᵏ_full = fₐᵏ_full / (α * ω₀²)
        fₓᵏ_full = fₓᵏ_full * β / (α^2 * ω₀²)
        Rₓ = Rₓ * β / (α^2 * ω₀²)
        Maxᵏ_full = Maxᵏ_full * (β / α)
        Mxaᵏ_full = Mxaᵏ_full * (β / α)
        ωₐ²ᵏ_full = ωₐ²ᵏ_full / ω₀²
        Mxxᵏ_full = Mxxᵏ_full * (β / α)^2
        Kxxᵏ_full = Kxxᵏ_full * (β / α)^2 * 1 / ω₀²



        kₙ, kₜ₁, kₜ₂ = kₙ * (β / α)^2 * 1 / ω₀², kₜ₁ * (β / α)^2 * 1 / ω₀²,
        kₜ₂ * (β / α)^2 * 1 / ω₀²

        γ = γ * β^4 / (α^2 * ω₀²)
        ξ = ξ * √(ω₀²)
    end

    fₐ = fₐᵏ_full[EO + 1]
    fₓ = fₓᵏ_full[EO + 1]

    force_factor = 1

    ext_force = ComplexHarmonicForcing([(1, ForcingVector(force_factor*fₐ, force_factor*fₓ))])
    # ext_force = HarmonicForcing([(1, ForcingVector(-im / 2 * 0.05fₐ, -im / 2 * 0.05fₓ))])

    cyclic_data = (
        ωₐ²ᵏ_full = ωₐ²ᵏ_full, Kxxᵏ_full = Kxxᵏ_full, Mxxᵏ_full = Mxxᵏ_full, Maxᵏ_full = Maxᵏ_full,
        Mxaᵏ_full = Mxaᵏ_full, fₐᵏ_full = fₐᵏ_full, fₓᵏ_full = fₓᵏ_full,
        Rₓ = Rₓ, Na = Na, Nx = Nx, Nc = Nc, Ns = Ns, Ne = Ne)


    cyclic_rom_struct = CyclicROMStructure(cyclic_data, EO, H, ext_force, ξ)


begin
    Ω₀, Ωₘₐₓ = 0.8, 1.2
    N = 2^7
    X_init = zeros((2H + 1) * Nx)
end


    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.02/EO, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-6, ϵₓ = 1e-6), verbose = true, ncols = 5)

    # friction_laws = [Coulomb_2_1D(kₙ = kₙ_vec[i], xₙ₀ = -0.268, kₜ₁ = kₜ₁_vec[i], kₜ₂ = kₜ₂_vec[i], μ = μ)
                    #    for i in 1:Nc]
      friction_laws = [Coulomb_2_1D(kₙ = kₙ, xₙ₀ = 0.0, kₜ₁ = kₜ₁, kₜ₂ = kₜ₂, μ = μ)
                       for i in 1:Nc]
    #  friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ
    Aⱼ = compute_cb_amplitudes_cyclic(sol)

    p = (
        EO, kₙ, kₜ₁, kₜ₂, μ, ξ, H, force_factor, ω₀², α, β, ContPars, N, extra_react,
        xp = hbm_prob.p.preload_state.xₚ
    )
    
    friction_name = ["1e-2", "25e-3", "5e-2", "1e-1", "15e-2"]
     @save "scripts/realistic_models/turbine_blade_model/study_of_friction/data/EO$(EO)/$(dataset[key])_mu$(friction_name[n]).jld2" Ω Xⱼ Aⱼ p

# end
# end

    A₁ = Aⱼ[1:(2H + 1), :]
    # A₁ = Aⱼ[(2H + 2):(4H + 2), :]
    # A₁ = Aⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^16, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    lines(Ω, a_max)


