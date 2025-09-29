using HBMContinuation, GLMakie, LinearAlgebra, TFG
use_formal_theme!()

    rom_data = read_rom_data_csv(joinpath(@__DIR__, "../../../data/rom/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    for i in eachindex(Rₓ)
        println(Rₓ[i])
    end

begin
    rom_data = read_rom_data_csv(joinpath(@__DIR__, "../../../data/rom/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    ξ = 1e-6
    kₙ, kₜ₁, kₜ₂ = 1e9, 1e6, 1e6
    μ=0.1

    ω₀² = ωₐ²[1]
    α = abs(fₐ[1] / ω₀²)
    β = abs(μ * Rₓ[end] / kₜ₁)

    fₐ = fₐ / (α * ω₀²)
    fₓ = fₓ * β / (α^2 * ω₀²)
    Rₓ = Rₓ * β / (α^2 * ω₀²)
    Max = Max * (β / α)
    ωₐ² = ωₐ² / ω₀²
    ξ = ξ * √(ω₀²)
    Mxx = Mxx * (β / α)^2
    Kxx = Kxx * (β / α)^2 * 1 / ω₀²
    kₙ, kₜ₁, kₜ₂ = kₙ * (β / α)^2 * 1 / ω₀², kₜ₁ * (β / α)^2 * 1 / ω₀²,
    kₜ₂ * (β / α)^2 * 1 / ω₀²

    ext_force = HarmonicForcing([(1, ForcingVector(fₐ, fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, Rₓ)
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Ω", ylabel = L"\text{max}|a_{1}|")
Hs = [1, 3, 5, 10]
use_formal_theme!()
# for i in eachindex(Hs)
begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 0.81, 0.85
    N, H = 2^7, Hs[i]
    X_init = zeros((2H + 1) * Nx)
end



    HBMCPars = HBMCParameters(N = N, H = H)

    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.01, maxsteps = 1_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-4, ϵₓ = 1e-6), verbose = true, ncols = 5)

function solver_HBM_coulomb(kₙ, xₙ₀, kₜ₁, kₜ₂, μ, rom_str, ContPars, HBMCPars, X_init, Ω₀)
    friction_laws = [Coulomb_2_1D(kₙ = kₙ, xₙ₀ = 0.0, kₜ₁ = kₜ₁, kₜ₂ = kₜ₂, μ = μ)
    for i in 1:Nc]

    hbm_prob = HBMCProblem(rom_str, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ

    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))
    return  a_max, Xⱼ, Ω 
end

    lines(Ω, a_max)
     #lines(Ω, a_max * α)

    X₁ = Xⱼ[1:(2H + 1), :]
     #X₁ = Xⱼ[(2H + 2):(4H + 2), :]
     #X₁ = Xⱼ[(4H + 3):(6H + 3), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    x₁ = E * X₁
    x_max = vec(maximum(abs.(x₁), dims = 1))

  # lines(Ω, x_max)
    # lines(Ω, x_max * β)

lines!(ax, Ω, a_max; label = "H = $(H)")

# end
axislegend(ax)
fig

using CairoMakie
save("scripts/tbm/1.tbm_simplificado/coulomb_diff_H.pdf", fig; backend=CairoMakie)
H = 5
using CairoMakie, Colors, ColorSchemes

function make_figure()
    use_formal_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")

    # ax.xticks = ([0.81, 0.82, 0.83, 0.84, 0.85], [L"$0.81$", L"$0.82$", L"$0.83$", L"$0.84$", L"$0.85$"])
    # ax.yticks = ([0, 50, 100, 150, 200, 250], [L"$0$", L"$50$", L"$100$", L"$150$", L"$200$", L"$250$"])

    return fig, ax
end

# Comparacion con diferentes fuerzas externas
begin
    factors = [0.5, 1.0, 2.0, 4.0]
    colors = [:blue, :red, :orange, :green]

    fig, ax = make_figure()

    for (i, factor) in enumerate(factors)
        ext_force = HarmonicForcing([(1, ForcingVector(factor * fₐ, factor * fₓ))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, Rₓ)

        a_max, Xⱼ, Ω = solver_HBM_coulomb(kₙ, 0.0, kₜ₁, kₜ₂, μ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"%$(factor) f_{\mathrm{ext}}")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures/fext_var_coulomb.pdf", fig)
end

#Comparacion con diferentes μs
begin
    μs = [0.15, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0]
    colors = distinguishable_colors(length(μs), ColorSchemes.tab10.colors)

    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")
    ax.xticks = ([0.75, 0.80, 0.85, 0.90], [L"$0.75$", L"$0.80$", L"$0.85$", L"$0.90$"])
    ax.yticks = ([0, 100, 200, 300], [L"$0$", L"$100$", L"$200$", L"$300$"])

    Ω₀, Ωₘₐₓ = 0.75, 0.90
    ContPars = ContinuationParameters(
        λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.01, maxsteps = 1_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(ϵᵣ = 1e-3, ϵₓ = 1e-5), verbose = true, ncols = 5)

    for (i, μᵢ) in enumerate(μs)
        ext_force = HarmonicForcing([(1, ForcingVector(fₐ, fₓ))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, Rₓ)

        a_max, Xⱼ, Ω = solver_HBM_coulomb(kₙ, 0.0, kₜ₁, kₜ₂, μᵢ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"μ = %$(μᵢ)")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/tbm/1.tbm_simplificado/mu_var_coulomb.pdf", fig)
end



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
lines(x₁[:, 350] .+ hbm_prob.p.preload_state.xₚ[3])

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
