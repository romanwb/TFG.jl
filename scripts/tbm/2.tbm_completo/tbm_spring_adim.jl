using  TFG, GLMakie, LinearAlgebra, Julianim, MAT
using HBMContinuation


    #  raw_data = read_rom_input("scripts/realistic_models/turbine_blade_model/data/Data_for_Dynamic_R0p75D.mat")
    #  reduced_data = reduce_rom(raw_data; n_modes=5)
    #  save_rom_data("scripts/realistic_models/turbine_blade_model/data/processed_data_2x2.mat", reduced_data)
    @unpack M, Max, Mxx, K, Kxx, ωₐ², fₐ, fₓ, Φ, Ψ, T_CB, Rₓ = load_reduced_rom("scripts/realistic_models/turbine_blade_model/data/processed_data_2x2.mat")
    ξ = 1e-6
    γ = 5e16
    Rₓ = zeros(length(fₓ))

    ω₀² = ωₐ²[1]
    α = abs(fₐ[1] / ω₀²)
    β = α

    γ = γ * β^4 / (α^2 * ω₀²)
    γ = 1e2
    fₐ = fₐ / (α * ω₀²)
    fₓ = fₓ * β / (α^2 * ω₀²)
    Rₓ = Rₓ * β / (α^2 * ω₀²)

    Max = Max * (β / α)
    ωₐ² = ωₐ² / ω₀²
    ξ = ξ * √(ω₀²)
    Mxx = Mxx * (β / α)^2
    Kxx = Kxx * (β / α)^2 * 1 / ω₀²
Kxx
    #! 
        Kxx = 0.5(Kxx+Kxx')
    #   Mxx = 0.5(Mxx+Mxx')
     Rₓ = vec(Rₓ)

    #! 
    
    ext_force = HarmonicForcing([(1, ForcingVector(1fₐ, 1fₓ))])
    rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

 println(Rₓ)

    Nx, Nc = rom_str.Nx, rom_str.Nc
    Ω₀, Ωₘₐₓ = 2300.0 / sqrt(ω₀²), 2500.0 / sqrt(ω₀²)
    N, H = 2^5, 2
    X_init = zeros((2H + 1) * Nx)



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
    # lines(Ω, a_max * α)


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
