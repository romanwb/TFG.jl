"Version original, usa datos ya procesados directamente del .csv"

using HBMContinuation, GLMakie, LinearAlgebra, Julianim

    # cyclic_rom_data = read_cyclic_data(joinpath(@__DIR__, "data/rom_cyclic/data_cyclic_CB.csv"))
    cyclic_rom_data = read_cyclic_rom(joinpath(@__DIR__, "data/rom_cyclic/data_cyclic_post.csv"))


    @unpack ωₐ²ᵏ_full, Kxxᵏ_full, Mxxᵏ_full, Maxᵏ_full,
    Mxaᵏ_full, fₐᵏ_full, fₓᵏ_full, Rₓ, Na, Nx, Nc, Ns, Ne = cyclic_rom_data

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
# writedlm("Omega_05.csv", Ω)
# writedlm("x_05.csv", a_max)

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
    # idx = 1
    # idx = 194
    # idx = 24
    idx = 40
    # idx = 278
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
