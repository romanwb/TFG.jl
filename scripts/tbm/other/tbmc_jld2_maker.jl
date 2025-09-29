using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
    # si aún no se ha convertido a .jld2:
    # convert_mat_to_jld2("data/tbm/ns146raw_eo0/Data_for_Dynamic_R1p5D.mat", "data/tbm/ns146raw_eo0/7x7_eo0.jld2")

    use_formal_theme!()
    @load "data/tbm/ns146_raw/ns146_first_last_rows.jld2" contenido

    begin
    Mee = contenido["Mee_new"]
    Mec = contenido["Mec_new"]
    Mcc = contenido["Mcc_new"]
    MCee = contenido["MCee_new"]
    Kee = contenido["Kee_new"]
    Kec = contenido["Kec_new"]
    Kcc = contenido["Kcc_new"]
    KCee = contenido["KCee_new"]
    f = sparse(contenido["F_vector"])
    fe = f[1:Int(sqrt(length(Mee))),:]
    fc = f[Int(sqrt(length(Mee)))+1:end,:]
    Rxe = spzeros(Int(sqrt(length(Mee))))
    Rₓ = sparsevec(contenido["Rc"])
    g = contenido["g"]
    xc = contenido["xc"]
    println(keys(contenido))
    norm(Kee)
    norm(Mee)
    norm(f)
    1e3
    Mee *= 1e3
    Mec *= 1e3
    Mcc *= 1e3
    MCee *= 1e3
    Kee *= 1e3
    Kec *= 1e3
    Kcc *= 1e3
    KCee *= 1e3
    end
    


     Rx = [Rxe; Rₓ]
    #  Rx = Rₓ
    # N=29
    # N = 29
    N = 146
    n_modes = 5

    ωₐ²ᵏ_full = Vector{Vector{Float64}}(undef, N)
    Maxᵏ_full = Vector{Matrix{ComplexF64}}(undef, N)
    Mxaᵏ_full = similar(Maxᵏ_full)
    Mxxᵏ_full = similar(Maxᵏ_full)
    Kxxᵏ_full = similar(Maxᵏ_full)
    fₐᵏ_full  = Vector{Vector{ComplexF64}}(undef, N)
    fₓᵏ_full  = similar(fₐᵏ_full)
    Φᵏ   = Vector{Matrix{ComplexF64}}(undef, N)
    Ψᵏ   = Vector{Matrix{ComplexF64}}(undef, N)

    EO = 1
    N = 146
    for k in 0:N-1
        _, Max, Mxx, _, Kxx, ω², fa, fx, Φ, Ψ, _ = craig_bampton_aplication_sym(Mee, Mec, Mcc, MCee, Kee, Kec, Kcc, KCee, fe, fc, Rx; k = k, N = N, EO=EO, n_modes = n_modes)
    
        Maxᵏ_full[k+1] = Max
        Mxaᵏ_full[k+1] = Max'
        Mxxᵏ_full[k+1] = Mxx
        Kxxᵏ_full[k+1] = Kxx
        fₐᵏ_full[k+1]  = vec(fa)
        fₓᵏ_full[k+1]  = vec(fx)
        ωₐ²ᵏ_full[k+1]  = ω²
        Φᵏ[k+1]   = Φ
        Ψᵏ[k+1]   = Ψ
    end

     @save "data/tbm/ns146_reduced/reduced_ns146_first_last_rows.jld2" Maxᵏ_full Mxaᵏ_full Mxxᵏ_full Kxxᵏ_full fₐᵏ_full fₓᵏ_full ωₐ²ᵏ_full Φᵏ Ψᵏ Rₓ


    #  @load "data/tbm/reduced_ns_29_EO0.jld2" Maxᵏ_full Mxaᵏ_full Mxxᵏ_full Kxxᵏ_full fₐᵏ_full fₓᵏ_full ωₐ²ᵏ_full
    # cyclic_rom_data = read_cyclic_rom(joinpath(@__DIR__, "data2/data_cyclic_CB.csv"))
    # @unpack Rₓ, Na, Nx, Nc, Ns, Ne = cyclic_rom_data
    #   Rₓ = [ -15.486772014026355, 2.449558123498594e-13, -50.41187612591457, 15.374436934430188, 2.55341482409175e-13, -50.30024088596551]
    Rₓ = Rx
    #   Rₓ = [ -5.2976896e+00, 8.3794083e-14, -1.7244812e+01, 5.2592622e+00, 8.7346796e-14, -1.7206624e+01]


    Na = 10
    Nx = 147
    Nc = 49
    Ns = 29
    # Ns = 8
    Ne = 36771

    kₙ, kₜ₁, kₜ₂, μ = 1e8, 2e6, 2e6, 0.1
    # kₙ, kₜ₁, kₜ₂, μ = 1e11, 8e10, 2e10, 0.1
    # kₙ, kₜ₁, kₜ₂, μ = 148871.069, 148871.069, 148871.069, 0.1
    ξ = 0.035
    γ = 5e16
    H = 5 # numero de armonicos

    ω₀² = ωₐ²ᵏ_full[EO + 1][1]
    α = abs(fₐᵏ_full[EO + 1][1] / ω₀²)
    # β = abs(0.1 * maximum(abs.(Rₓ)) / kₜ₁)
    β = α

    #  ξ_vals = [0.0001, 0.035, 1.0]
    # begin
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
    # end
    fₐ = fₐᵏ_full[EO + 1]
    fₓ = fₓᵏ_full[EO + 1]



    ext_force = ComplexHarmonicForcing([(1, ForcingVector(fₐ, fₓ))])
    # ext_force = HarmonicForcing([(1, ForcingVector(-im / 2 * 0.05fₐ, -im / 2 * 0.05fₓ))])

    cyclic_data = (
        ωₐ²ᵏ_full = ωₐ²ᵏ_full, Kxxᵏ_full = Kxxᵏ_full, Mxxᵏ_full = Mxxᵏ_full, Maxᵏ_full = Maxᵏ_full,
        Mxaᵏ_full = Mxaᵏ_full, fₐᵏ_full = fₐᵏ_full, fₓᵏ_full = fₓᵏ_full,
        Rₓ = Rₓ, Na = Na, Nx = Nx, Nc = Nc, Ns = Ns, Ne = Ne)


    cyclic_rom_struct = CyclicROMStructure(cyclic_data, EO, H, ext_force, ξ)



begin
    Ω₀, Ωₘₐₓ = 0.8, 1.2
    N = 2^9
    X_init = zeros((2H + 1) * Nx)
end

begin
    # for ξ in ξ_vals
    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.005, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-6, ϵₓ = 1e-6), verbose = true, ncols = 5)

    # friction_laws = [Coulomb_2_1D(kₙ = kₙ_vec[i], xₙ₀ = -0.268, kₜ₁ = kₜ₁_vec[i], kₜ₂ = kₜ₂_vec[i], μ = μ)
                    #    for i in 1:Nc]
    #   friction_laws = [Coulomb_2_1D(kₙ = kₙ, xₙ₀ = 0.0, kₜ₁ = kₜ₁, kₜ₂ = kₜ₂, μ = μ)
                    #    for i in 1:Nc]
     friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws)
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

    # lines(Ω, x_max)

    A₁ = Aⱼ[1:(2H + 1), :]
    # A₁ = Aⱼ[(2H + 2):(4H + 2), :]
    # A₁ = Aⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    lines(Ω, a_max)
    # # lines(Ω, a_max * α)
    # end
end

# ================================#
#=================================#


begin
    idx = 1450
    c_idx = 20

    X₁ = Xⱼ[(1 + (2H + 1) * 3(c_idx - 1)):(1 + 2H + (2H + 1) * 3(c_idx - 1)), :]
    X₂ = Xⱼ[(2H + 2 + (2H + 1) * 3(c_idx - 1)):(4H + 2 + (2H + 1) * 3(c_idx - 1)), :]
    X₃ = Xⱼ[(4H + 3 + (2H + 1) * 3(c_idx - 1)):(6H + 3 + (2H + 1) * 3(c_idx - 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    xp = hbm_prob.p.preload_state.xₚ
    x₁ = (E * X₁)[:, idx] .+ xp[1 + 3 * (c_idx - 1)]
    x₂ = (E * X₂)[:, idx] .+ xp[2 + 3 * (c_idx - 1)]
    x₃ = (E * X₃)[:, idx] .+ xp[3 + 3 * (c_idx - 1)]

    t_force1 = zeros(2^10)
    t_force2 = zeros(2^10)
    normal = HBMContinuation.normal_force.(x₃, kₙ, 0.0)
    ws1, ws2 = 0.0, 0.0
    for j in 1:2
        for i in eachindex(t_force1)
            t_force2[i], w2 = HBMContinuation.tangential_force(
                x₂[i], ws2, kₜ₂, μ, normal[i])
            ws2 = w2

            t_force1[i], w1 = HBMContinuation.tangential_force(
                x₁[i], ws1, kₜ₁, μ, normal[i])
            ws1 = w1
        end
    end

    lines(x₁, t_force1)
      lines(x₂, t_force2)
end


    # idx = 1
    # idx = 194
    idx = 1450
    c_idx = 1
    X₁ = Xⱼ[(1 + (2H + 1) * 3(c_idx - 1)):(1 + 2H + (2H + 1) * 3(c_idx - 1)), :]
    X₂ = Xⱼ[(2H + 2 + (2H + 1) * 3(c_idx - 1)):(4H + 2 + (2H + 1) * 3(c_idx - 1)), :]
    X₃ = Xⱼ[(4H + 3 + (2H + 1) * 3(c_idx - 1)):(6H + 3 + (2H + 1) * 3(c_idx - 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    xp = hbm_prob.p.preload_state.xₚ
    x₁ = (E * X₁)[:, idx] .+ xp[1 + 3 * (c_idx - 1)]
    x₂ = (E * X₂)[:, idx] .+ xp[2 + 3 * (c_idx - 1)]
    x₃ = (E * X₃)[:, idx] .+ xp[3 + 3 * (c_idx - 1)]

    n_cycles = 2

    t_force1 = zeros(2^10 * n_cycles)
    t_force2 = zeros(2^10 * n_cycles)
    normal = HBMContinuation.normal_force.(x₃, kₙ, 0.0)
    ws1, ws2 = 0.0, 0.0
    for j in 1:n_cycles
        for i in eachindex(x₁)
            t_force2[i + 2^10 * (j - 1)], w2 = HBMContinuation.tangential_force(
                x₂[i], ws2, kₜ₂, μ, normal[i])
            ws2 = w2

            t_force1[i + 2^10 * (j - 1)], w1 = HBMContinuation.tangential_force(
                x₁[i], ws1, kₜ₁, μ, normal[i])
            ws1 = w1
        end
    end

    t = range(0, stop = 2π * n_cycles, length = 2^10 * n_cycles)

    x₁ = repeat(x₁, n_cycles)
    x₂ = repeat(x₂, n_cycles)

    lines(x₁, t_force1)
    lines(x₂, t_force2)

    lines(t, x₁)
    lines(t, t_force1)
    lines(t, x₂)
    lines(t, t_force2)


begin
    #  idx = 1450
    #  idx = 194
    # idx = 24
    idx = 40
    #  idx = 278
    c_idx = 11
    X₁ = Xⱼ[(1 + (2H + 1) * 3(c_idx - 1)):(1 + 2H + (2H + 1) * 3(c_idx - 1)), :]
    X₂ = Xⱼ[(2H + 2 + (2H + 1) * 3(c_idx - 1)):(4H + 2 + (2H + 1) * 3(c_idx - 1)), :]
    X₃ = Xⱼ[(4H + 3 + (2H + 1) * 3(c_idx - 1)):(6H + 3 + (2H + 1) * 3(c_idx - 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    xp = hbm_prob.p.preload_state.xₚ
    x₁ = (E * X₁)[:, idx] .+ xp[1 + 3 * (c_idx - 1)]
    x₂ = (E * X₂)[:, idx] .+ xp[2 + 3 * (c_idx - 1)]
    x₃ = (E * X₃)[:, idx] .+ xp[3 + 3 * (c_idx - 1)]

    n_cycles = 5

    t_force1 = zeros(2^10 * n_cycles)
    t_force2 = zeros(2^10 * n_cycles)
    normal = HBMContinuation.normal_force.(x₃, kₙ, 0.0)
    w1_vec = zeros(2^10 * n_cycles)
    w2_vec = zeros(2^10 * n_cycles)
    ws1, ws2 = 0.0, 0.0
    for j in 1:n_cycles
        for i in eachindex(x₁)
            t_force2[i + 2^10 * (j - 1)], w2 = HBMContinuation.tangential_force(
                x₂[i], ws2, kₜ₂, μ, normal[i])
            ws2 = w2
            w2_vec[i + 2^10 * (j - 1)] = w2

            t_force1[i + 2^10 * (j - 1)], w1 = HBMContinuation.tangential_force(
                x₁[i], ws1, kₜ₁, μ, normal[i])
            ws1 = w1
            w1_vec[i + 2^10 * (j - 1)] = w1
        end
    end

    t = range(0, stop = 2π * n_cycles, length = 2^10 * n_cycles)

    x₁ = repeat(x₁, n_cycles)
    x₂ = repeat(x₂, n_cycles)
end
    lines(x₁, t_force1)
    lines(x₂, t_force2)

    # lines(t, x₁)
    # lines(t, t_force1)
    # lines(t, x₂)
    # lines(t, t_force2)
    # lines(t, w1_vec)
    # lines(t, w2_vec)

