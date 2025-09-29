using HBMContinuation, GLMakie, LinearAlgebra, Julianim

begin
    rom_data = read_rom_data_csv(joinpath(@__DIR__, "../../../data/rom/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    ξ = 1e-6

    ext_force = HarmonicForcing([(1, ForcingVector(fₐ, fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, Rₓ)
end


begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 3800.0, 5000.0
    N, H = 2^7, 10
    X_init = zeros((2H + 1) * Nx)
end

begin
    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 1.0, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-4, ϵₓ = 1e-6), verbose = true, ncols = 5)

    friction_laws = [Coulomb_2_1D(kₙ = 1e9, xₙ₀ = 0.0, kₜ₁ = 1e6, kₜ₂ = 2e6, μ = 0.01)
                     for i in 1:Nc]

    hbm_prob = HBMCProblem(rom_str, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ

    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))
    lines(Ω, a_max)

    X₁ = Xⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))
end

lines(Ω, a_max)

using DelimitedFiles
writedlm("Omega10.csv", Ω)
writedlm("x10.csv", a_max)

begin
    X₁ = Xⱼ[(0(2H + 1) + 1):(1(2H + 1)), :]
    # X₁ = Xⱼ[(2(2H + 1) + 1):(3(2H + 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))
    lines(Ω, x_max)
end

lines(x₁[:, 200] .+ hbm_prob.p.preload_state.xₚ[1])
# lines(x₁[:, 350] .+ hbm_prob.p.preload_state.xₚ[3])

# Ω = 4977 229, Ω = 5000, 321

# using DelimitedFiles
# M = [I Max;
#      Max' Mxx]

# K = [diagm(ωₐ²) zeros(5, Nx);
#      zeros(Nx, 5) Kxx]

# f = [fₐ; fₓ]

# xp = hbm_prob.p.preload_state.xₚ
# gp = hbm_prob.p.preload_state.gₚ

# writedlm("M.csv", M)
# writedlm("K.csv", K)
# writedlm("f.csv", f)
# writedlm("xp.csv", xp)
# writedlm("gp.csv", gp)

begin
    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1])

    ax.xticks = ([0, π / 2, π, 3π / 2, 2π],
        [L"$0$", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"])

    t = range(0, 2π, length = length(a₁[:, 476]))

    lines!(ax, [0, 2π], [0, 0]; color = :black, linestyle = :dash, linewidth = 1)
    lines!(ax, t, x₁[:, 321]; color = :black, linewidth = 6)

    ax.xlabel = L"$\Omega t$"
    ax.ylabel = L"$x_1(t)$"

    # ylims!(ax, (-0.55, 0.55))

    # save("x_1_5000.png", fig)

    display(GLMakie.Screen(), fig)
end
