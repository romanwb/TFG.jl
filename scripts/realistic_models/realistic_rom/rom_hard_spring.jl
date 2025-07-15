using HBMContinuation, GLMakie, LinearAlgebra, Julianim

begin
    rom_data = read_rom_data_csv(joinpath(@__DIR__, "data/data.csv"))
    @unpack Max, Mxx, Kxx, ωₐ², fₐ, fₓ, Rₓ = rom_data
    ξ = 1e-6

    ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)
end

begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 4920.0, 5050.0
    γ = 5e16
    N, H = 2^5, 10
    X_init = zeros((2H + 1) * Nx)
end


HBMCPars = HBMCParameters(N = N, H = H)

ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.25, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

function solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)
    friction_laws = [Cubic(γ, γ, γ) for i in 1:Nc]

    hbm_prob = HBMCProblem(rom_str, ContPars, HBMCPars, friction_laws; finite_diff = false)
    sol = continuation(hbm_prob, X_init, Ω₀)

    Xⱼ, Ω = sol.x, sol.λ

    Aⱼ = compute_cb_amplitudes(sol)

    A₁ = Aⱼ[1:(2H + 1), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    a₁ = E * A₁
    a_max = vec(maximum(abs.(a₁), dims = 1))

    return a_max, Xⱼ, Ω
end

a_max, Xⱼ, Ω = solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

X₁ = Xⱼ[(0(2H + 1) + 1):(1(2H + 1)), :]
E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
x₁ = E * X₁

using CairoMakie

function make_figure()
    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")

    ax.xticks = ([4950, 5000, 5050], [L"$4950$", L"$5000$", L"$5050$"])
    ax.yticks = ([2.5e-5, 5e-5, 7.5e-5], [L"$2.5 \times 10^{-5}$", L"$5.0 \times 10^{-5}$", L"$7.5 \times 10^{-5}$"])

    return fig, ax
end



# Comparacion con diferente ξ
begin
    ξ_vals = [1.0, 1.25, 1.5, 1.75, 2.0]
    colors = [:blue, :red, :orange, :green, :purple]

    fig, ax = make_figure()

    for (i, ξ) in enumerate(ξ_vals)
        ξ
        ξ_real = ξ * 1e-6
        ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ_real, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"\xi = %$(ξ) \text{e} -6")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures/xi_var.pdf", fig)
end


# Comparacion con diferentes fuerzas externas
begin
    factors = [1.0, 0.75, 0.5, 0.25]
    colors = [:blue, :red, :orange, :green]

    fig, ax = make_figure()

    for (i, factor) in enumerate(factors)
        ext_force = HarmonicForcing([(1, ForcingVector(factor * fₐ, factor * fₓ))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"%$(factor) f_{\mathrm{ext}}")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures/fext_var.pdf", fig)
end

# Comparacion con diferente γ
begin
    γ_vals = [5e16, 5e17, 5e18, 5e20]
    colors = [:blue, :red, :orange, :green]
    Ω₀, Ωₘₐₓ = 4900.0, 5200.0

    ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.25, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], yscale = log10, xlabel = L"\Omega", ylabel = L"|a_1(t)|")
    ax.xticks = ([4900, 5025, 5150], [L"$4900$", L"$5025$", L"$5150$"])
    ax.yticks = ([10^(-5.5), 10^(-5), 10^(-4.5), 10^(-4)], [L"10^{-5.5}", L"10^{-5}", L"10^{-4.5}", L"10^{-4.5}"])

    for (i, γᵢ) in enumerate(γ_vals)
        ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])

        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γᵢ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        label_str = L"\gamma = %$(γᵢ)"
        lines!(ax, Ω, a_max; color = colors[i], label = label_str)
    end

    vlines!(ax, [5078.0]; color = :black, linestyle = :dashdot, linewidth = 1)

    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures/gamma_var.pdf", fig)
end


