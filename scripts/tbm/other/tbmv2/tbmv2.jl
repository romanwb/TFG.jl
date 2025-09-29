using  TFG, GLMakie, LinearAlgebra, Julianim, MAT, HBMContinuation

@unpack M, Max, Mxx, K, Kxx, ωₐ², fₐ, fₓ, Φ, Ψ, T_CB = load_reduced_rom("data/tbm/processed_data_7x7.mat")

ξ = 1e-6
ext_force = HarmonicForcing([(1, ForcingVector(1*real(fₐ), 1*real(fₓ)))])
Rₓ = ones(length(fₓ))
rom_str = ROMStructure(Max, Mxx, Kxx, ωₐ², ext_force, ξ, 10Rₓ)

sqrt(ωₐ²[1])
begin
    Nx, Nc = rom_str.Nx, rom_str.Nc
    #Ω₀, Ωₘₐₓ = 4920.0, 5050.0
    Ω₀, Ωₘₐₓ = 2300.0, 2500.0
    γ = 5e20
    N, H = 2^5, 10
    X_init = zeros((2H + 1) * Nx)
end

HBMCPars = HBMCParameters(N = N, H = H)

ContPars = ContinuationParameters(λmin = Ω₀, λmax = Ωₘₐₓ, Δs = 5.0, maxsteps = 10_000,
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
