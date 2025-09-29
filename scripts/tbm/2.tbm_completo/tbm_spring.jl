using  TFG, GLMakie, LinearAlgebra, Julianim, MAT, HBMContinuation, JLD2

@load "data/tbm/ns146raw_eo0/7x7_eo0_processed.jld2" Max Mxx Kxx fₐ fₓ ω² Φ Ψ Rx
Mxa = Max'
Nx = 147

Rₓ = Rx[end-146:end]

ξ = 1e-6
kₙ, kₜ₁, kₜ₂, μ = 1e9, 1e6, 1e6, 0.01
γ = 5e16
#adim
begin
    # ω₀² = 2.0456228983408358e7
    ω₀² = real(ω²[1])
    α = abs(fₐ[1] / ω₀²)
    # β = abs(0.1 * maximum(abs.(Rₓ)) / kₜ₁)
    β = α
    fₐ = fₐ / (α * ω₀²)
    fₓ = fₓ * β / (α^2 * ω₀²)
    Rₓ = Rₓ * β / (α^2 * ω₀²)
    Max = Max * (β / α)
    Mxa = Mxa * (β / α)
    ω² = ω² / ω₀²
    Mxx = Mxx * (β / α)^2
    Kxx = Kxx * (β / α)^2 * 1 / ω₀²
    kₙ, kₜ₁, kₜ₂ = kₙ * (β / α)^2 * 1 / ω₀², kₜ₁ * (β / α)^2 * 1 / ω₀²,
    kₜ₂ * (β / α)^2 * 1 / ω₀²
    γ = γ * β^4 / (α^2 * ω₀²)
    ξ = ξ * √(ω₀²)
end

ext_force = HarmonicForcing([(1, ForcingVector(1*real(vec(fₐ)), 1*real(vec(fₓ))))])
# ext_force = HarmonicForcing([(1, ForcingVector(1*real(vec(fₐ)), 1*real(vec(fₓ)))), (2, ForcingVector(50*real(vec(fₐ)), 50*real(vec(fₓ))))])
rom_str = ROMStructure(Max, Mxx, Kxx, ω², ext_force, ξ, Rₓ)

sqrt(ω²[1])
begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    #Ω₀, Ωₘₐₓ = 4920.0, 5050.0
    # Ω₀, Ωₘₐₓ = 2300.0, 2500.0
    Ω₀, Ωₘₐₓ = 0.97, 1.03
    N, H = 2^7, 1
    X_init = zeros((2H + 1) * Nx)
end

HBMCPars = HBMCParameters(N = N, H = H)

ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.001, maxsteps = 10_000,
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

lines(Ω, a_max)


# using CairoMakie

function make_figure()
    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")

    ax.xticks = ([2300, 2400, 2500], [L"$2300$", L"$2400$", L"$2500$"])
    ax.yticks = ([1e-7, 2e-7, 3e-7], [L"$1 \times 10^{-7}$", L"$2 \times 10^{-7}$", L"$3 \times 10^{-7}$"])

    return fig, ax
end


# Comparacion con diferente ξ
begin
    ξ_vals = [1.0, 1.25, 1.5, 1.75, 2.0]
    colors = [:blue, :red, :orange, :green, :purple]

    # fig, ax = make_figure()

    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")

    for (i, ξ) in enumerate(ξ_vals)
        ξ
        ξ_real = ξ * 1e-6
        ext_force = HarmonicForcing([(1, ForcingVector(1*real(fₐ), 1*real(fₓ)))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ_real, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"\xi = %$(ξ) \text{e} -6")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures_matlab_problem/xi_var.pdf", fig)
end


# Comparacion con diferentes fuerzas externas
begin
    factors = [1.0, 0.75, 0.5, 0.25]
    colors = [:blue, :red, :orange, :green]

    # fig, ax = make_figure()

    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"|a_1(t)|")

    for (i, factor) in enumerate(factors)
        ext_force = HarmonicForcing([(1, ForcingVector(factor * real(fₐ), factor * real(fₓ)))])
        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        lines!(ax, Ω, a_max; color = colors[i], label = L"%$(factor) f_{\mathrm{ext}}")
    end
    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures_matlab_problem/fext_var.pdf", fig)
end

# Comparacion con diferente γ
begin
    # γ_vals = [5e16, 5e17, 5e18, 5e20]
    γ_vals = [5e1, 5e3, 5e10, 5e20]
    colors = [:blue, :red, :orange, :green]

    ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 0.25, maxsteps = 10_000,
        direction = :forward, predictor = PseudoArcLength(),
        corrector = Newton(), verbose = true, ncols = 5)

    set_publication_theme!()
    fig = Figure()
    ax = Axis(fig[1, 1], yscale = log10, xlabel = L"\Omega", ylabel = L"|a_1(t)|")
    # ax.xticks = ([4900, 5025, 5150], [L"$4900$", L"$5025$", L"$5150$"])
    # ax.yticks = ([10^(-5.5), 10^(-5), 10^(-4.5), 10^(-4)], [L"10^{-5.5}", L"10^{-5}", L"10^{-4.5}", L"10^{-4.5}"])

    for (i, γᵢ) in enumerate(γ_vals)
        ext_force = HarmonicForcing([(1, ForcingVector(1*real(fₐ), 1*real(fₓ)))])

        rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

        a_max, Xⱼ, Ω = solver_HBM(γᵢ, rom_str, ContPars, HBMCPars, X_init, Ω₀)

        label_str = L"\gamma = %$(γᵢ)"
        lines!(ax, Ω, a_max; color = colors[i], label = label_str)
    end

    # vlines!(ax, [5078.0]; color = :black, linestyle = :dashdot, linewidth = 1)

    axislegend(ax; position = :rt)

    display(GLMakie.Screen(), fig)
    save("scripts/rom_model/figures_matlab_problem/gamma_var.pdf", fig)
end
