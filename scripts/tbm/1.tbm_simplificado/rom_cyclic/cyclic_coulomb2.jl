"Versión completa, que usa los datos en crudo del .txt (previamente convertidos a .jld2), hace transformación de CB,..."

using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2

    # cyclic_rom_data = read_cyclic_data(joinpath(@__DIR__, "data/data_cyclic_CB.csv"))
    # cyclic_rom_data = read_cyclic_rom(joinpath(@__DIR__, "data2/data_cyclic_CB.csv"))

    # Se cargan los datos guardados en el JLD2
    @load "data/rom_cyclic/FEM_matrices24-6.jld2" datos
    Mee = datos["Mee"]
    Mec = datos["Mec"]
    Mcc = datos["Mcc"]
    MCee = datos["MCee"]
    Kee = datos["Kee"]
    Kec = datos["Kec"]
    Kcc = datos["Kcc"]
    KCee = datos["KCee"]
    fe = datos["Fe"]
    fc = datos["Fc"]
    Rxe = datos["Pe"]
    Rxc = datos["Pc"]
    Rx = [Rxe, Rxc]

    N=8
    n_modes = 10

    ωₐ²ᵏ_full = Vector{Vector{Float64}}(undef, N)

    Maxᵏ_full = Vector{Matrix{ComplexF64}}(undef, N)
    Mxaᵏ_full = similar(Maxᵏ_full)
    Mxxᵏ_full = similar(Maxᵏ_full)
    Kxxᵏ_full = similar(Maxᵏ_full)
    fₐᵏ_full  = Vector{Vector{ComplexF64}}(undef, N)
    fₓᵏ_full  = similar(fₐᵏ_full)
    # Φᵏ   = Vector{Matrix{ComplexF64}}(undef, N)
    # Ψᵏ   = Vector{Matrix{ComplexF64}}(undef, N)

    for k in 0:N-1
        _, Max, Mxx, _, Kxx, ω², fa, fx, Φ, Ψ, _ = craig_bampton_aplication_sym(Mee, Mec, Mcc, MCee, Kee, Kec, Kcc, KCee, fe, fc, Rx; k = k, N = N, EO=1, n_modes = n_modes)
    
        Maxᵏ_full[k+1] = Max
        Mxaᵏ_full[k+1] = Max'
        Mxxᵏ_full[k+1] = Mxx
        Kxxᵏ_full[k+1] = Kxx
        fₐᵏ_full[k+1]  = vec(fa)
        fₓᵏ_full[k+1]  = vec(fx)
        ωₐ²ᵏ_full[k+1]  = ω²
        # Φᵏ[k+1]   = Φ
        # Ψᵏ[k+1]   = Ψ
    end

    # cyclic_rom_data = read_cyclic_rom(joinpath(@__DIR__, "data2/data_cyclic_CB.csv"))
    # @unpack Rₓ, Na, Nx, Nc, Ns, Ne = cyclic_rom_data
    #   Rₓ = [ -15.486772014026355, 2.449558123498594e-13, -50.41187612591457, 15.374436934430188, 2.55341482409175e-13, -50.30024088596551]

     Rₓ = [ -5.2976896e+00, 8.3794083e-14, -1.7244812e+01, 5.2592622e+00, 8.7346796e-14, -1.7206624e+01]


    Na = 10
    Nx = 6
    Nc = 2
    Ns = 8
    Ne = 180

    kₙ, kₜ₁, kₜ₂, μ = 1e8, 2e6, 2e6, 0.1

    β = abs(0.1 * maximum(abs.(Rₓ)) / kₜ₁)

    H = 5
    EO = 1

    ω₀² = ωₐ²ᵏ_full[EO + 1][1]
    α = abs(fₐᵏ_full[EO + 1][1] / ω₀²)

    ξ = 0.035

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

    fₐ = fₐᵏ_full[EO + 1]
    fₓ = fₓᵏ_full[EO + 1]

    ext_force = ComplexHarmonicForcing([(1, ForcingVector(0.01fₐ, 0.01fₓ))])
    # ext_force = HarmonicForcing([(1, ForcingVector(-im / 2 * 0.05fₐ, -im / 2 * 0.05fₓ))])

    cyclic_data = (
        ωₐ²ᵏ_full = ωₐ²ᵏ_full, Kxxᵏ_full = Kxxᵏ_full, Mxxᵏ_full = Mxxᵏ_full, Maxᵏ_full = Maxᵏ_full,
        Mxaᵏ_full = Mxaᵏ_full, fₐᵏ_full = fₐᵏ_full, fₓᵏ_full = fₓᵏ_full,
        Rₓ = Rₓ, Na = Na, Nx = Nx, Nc = Nc, Ns = Ns, Ne = Ne)

    cyclic_rom_struct = CyclicROMStructure(cyclic_data, EO, H, ext_force, ξ)


begin
    Ω₀, Ωₘₐₓ = 0.35, 0.45
    N = 2^9
    X_init = zeros((2H + 1) * Nx)
end

begin
    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.005, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-8, ϵₓ = 1e-8), verbose = true, ncols = 5)

    friction_laws = [Coulomb_2_1D(kₙ = kₙ, xₙ₀ = 0.0, kₜ₁ = kₜ₁, kₜ₂ = kₜ₂, μ = μ)
                     for i in 1:Nc]

    hbm_prob = CyclicHBMCProblem(cyclic_rom_struct, ContPars, HBMCPars, friction_laws)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ
    Aⱼ = compute_cb_amplitudes_cyclic(sol)

    X₁ = Xⱼ[1:(2H + 1), :]
    # X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))

    lines(Ω, x_max)

    # A₁ = Aⱼ[1:(2H + 1), :]
    # # A₁ = Aⱼ[(2H + 2):(4H + 2), :]
    # # A₁ = Aⱼ[(4H + 3):(6H + 3), :]
    # E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    # a₁ = E * A₁
    # a_max = vec(maximum(abs.(a₁), dims = 1))

    # lines(Ω, a_max)
    # # lines(Ω, a_max * α)
end

# using DelimitedFiles
# writedlm("Omega10.csv", Ω)
# writedlm("x10.csv", a_max)

begin
    Xⱼ, Ω = sol.x, sol.λ

    X₁ = Xⱼ[1:(2H + 1), :]
    # X₁ = Xⱼ[(2H + 2):(4H + 2), :]
    # X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))

    lines(Ω, x_max)
    # lines(Ω, x_max * β)
end

begin
    idx = 170
    c_idx = 1

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
    # lines(x₂, t_force2)
end

# using DelimitedFiles
# # writedlm("Omega_05.csv", Ω)
# # writedlm("x_05.csv", a_max)

# writedlm("x_mu0025.csv", x₂)
# writedlm("T_mu0025.csv", t_force2)

begin
    # idx = 1
    # idx = 194
    idx = 40
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
end

# writedlm("t_5cycles.csv", t)

# writedlm("x1_5cycles.csv", x₁)
# writedlm("T1_5cycles.csv", t_force1)

# writedlm("x2_5cycles.csv", x₂)
# writedlm("T2_5cycles.csv", t_force2)

begin
     idx = 1
    #  idx = 194
    # idx = 24
    # idx = 40
    #  idx = 278
    c_idx = 2
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

    lines(x₁, t_force1)
    lines(x₂, t_force2)

    lines(t, x₁)
    lines(t, t_force1)
    lines(t, x₂)
    lines(t, t_force2)
    lines(t, w1_vec)
    lines(t, w2_vec)
end

# writedlm("t_5cycles.csv", t)

# writedlm("x1_5cycles.csv", x₁)
# writedlm("w1_5cycles.csv", w1_vec)
# writedlm("T1_5cycles.csv", t_force1)

# writedlm("x2_5cycles.csv", x₂)
# writedlm("w2_5cycles.csv", w2_vec)
# writedlm("T2_5cycles.csv", t_force2)
