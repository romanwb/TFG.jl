using HBMContinuation, LinearAlgebra, Julianim

begin
    rom_data = read_rom_data_csv(joinpath(@__DIR__, "data/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    ξ = 1e-6
    γ = 5e16

    ω₀² = ωₐ²[1]
    α = abs(fₐ[1] / ω₀²)
    β = α

    γ = γ * β^4 / (α^2 * ω₀²)
    fₐ = fₐ / (α * ω₀²)
    fₓ = fₓ * β / (α^2 * ω₀²)
    Rₓ = Rₓ * β / (α^2 * ω₀²)
    Max = Max * (β / α)
    ωₐ² = ωₐ² / ω₀²
    ξ = ξ * √(ω₀²)
    Mxx = Mxx * (β / α)^2
    Kxx = Kxx * (β / α)^2 * 1 / ω₀²

    ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)
end
Kxx
begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 4920.0 / sqrt(ω₀²), 5050.0 / sqrt(ω₀²)
    N, H = 2^5, 2
    X_init = zeros((2H + 1) * Nx)
end

begin
    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.1, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

    friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = HBMCProblem(rom_str, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ

    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    lines(Ω, a_max)
    # lines(Ω, a_max * α)
end

X₁ = Xⱼ[(0(2H + 1) + 1):(1(2H + 1)), :]
E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
x₁ = E * X₁

# Ω = 4977 229, Ω = 5000, 321

using DelimitedFiles
M = [I Max;
     Max' Mxx]

K = [diagm(ωₐ²) zeros(5, Nx);
     zeros(Nx, 5) Kxx]

f = [fₐ; fₓ]

xp = hbm_prob.p.preload_state.xₚ
gp = hbm_prob.p.preload_state.gₚ

writedlm("M.csv", M)
writedlm("K.csv", K)
writedlm("f.csv", f)
writedlm("xp.csv", xp)
writedlm("gp.csv", gp)

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
