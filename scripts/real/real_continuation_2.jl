# ========================== FUNCIONES AUXILIARES =============================
using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

function assemble_M_and_K(n::Int, m::Float64, k::Float64)
    M = Diagonal(fill(m, n)) |> Matrix
    K = zeros(n, n)
    for i in 1:(n - 1)
        K[i, i] += k
        K[i+1, i+1] += k
        K[i, i+1] -= k
        K[i+1, i] -= k
    end
    return M, K
end

assemble_C(_, K, ξ) = ξ * K

function define_forcing_vector(n::Int, Fp::Vector, t::Vector)
    @assert length(Fp) == n "Length of precarga Fp must match number of DOFs"
    Nt = length(t)
    F = zeros(n, Nt)
    for i in 1:n
        F[i, :] .= Fp[i] .* (sin.(t) .+ sin.(2 .* t))
    end
    return F
end

function apply_fft_projection(Eᴴ, F)
    return Eᴴ * F'
end

function apply_fourier_matrix(E, x̂::Vector, n::Int, Nh::Int)
    reshape(vcat([E * x̂[(i-1)*Nh+1:i*Nh] for i in 1:n]...), n * size(E, 1))
end

function apply_fourier_adjoint(Eᴴ, X::Vector, n::Int, Nh::Int)
    reshape(vcat([Eᴴ * X[(i-1)*size(Eᴴ, 2)+1:i*size(Eᴴ, 2)] for i in 1:n]...), n * Nh)
end

function define_nonlinearity(n::Int, nonlinear_dofs::Vector{Int}; type=:cubic)
    g(x) = begin
        T = eltype(x)
        out = zeros(T, n)
        if type == :cubic
            out[nonlinear_dofs] .= x[nonlinear_dofs] .^ 3
        end
        return out
    end
    dg(x) = begin
        T = eltype(x)
        out = zeros(T, n)
        if type == :cubic
            out[nonlinear_dofs] .= 3 .* x[nonlinear_dofs] .^ 2
        end
        return out
    end
    return g, dg
end

function build_A_hbm(M, C, K, H, ω)
    n = size(M, 1)
    Nh = 2H + 1
    T = typeof(ω)
    A = zeros(T, n * Nh, n * Nh)

    A[1:n, 1:n] .= Matrix{T}(I, n, n)

    for j in 1:H
        jω = T(j) * ω
        idx_b = n * (2j - 1) + 1
        idx_a = n * (2j) + 1
        Ab = - (jω)^2 * M + jω * C + K
        Aa = - (jω)^2 * M - jω * C + K
        A[idx_b:idx_b+n-1, idx_b:idx_b+n-1] .= Ab
        A[idx_a:idx_a+n-1, idx_a:idx_a+n-1] .= Aa
    end
    return A
end


function equilibrium_to_fourier(x0::Vector, n::Int, H::Int)
    Nh = 2H + 1
    x̂₀ = zeros(n * Nh)
    x̂₀[1:n] .= x0
    return x̂₀
end

function residual_hbm(x̂, λ, p)
    E, Eᴴ, H = p.E, p.Eᴴ, p.H
    M, C, K = p.M, p.C, p.K
    f̂, g = p.f̂, p.g

    n = size(M, 1)
    Nh = 2H + 1
    T = eltype(x̂)

    Xreal = apply_fourier_matrix(E, x̂, n, Nh)
    Xreal_mat = reshape(Xreal, n, :)
    gX_flat = reduce(vcat, [g(Xreal_mat[:, j]) for j in 1:size(Xreal_mat, 2)])
    ĝ = apply_fourier_adjoint(Eᴴ, gX_flat, n, Nh)

    A = build_A_hbm(M, C, K, H, λ)
    res = A * x̂ + ĝ - f̂
    return res
end

# ========================== CONFIGURACIÓN Y RESOLUCIÓN ======================
N, H = 2^6, 3
n = 5
Nh = 2H + 1
ξ, ϵ = 0.05, 1.0
ξ̃ = ϵ * ξ
λ₀ = 0.0

E, Eᴴ = fft_matrices(N, H)
M, K = assemble_M_and_K(n, 1.0, 10.0)
K[1, 1] += 1e-3  # penalización suave para evitar singularidad
C = assemble_C(M, K, ξ)

Fp_static = [1.0, 0.1, 0.1, 0.1, 0.1]
Fp_periodic = [1.0, 1.0, 0.0, 0.0, 0.0]

t = collect(range(0, 2π, length = N + 1)[1:end-1])
F_static = repeat(Fp_static, 1, length(t))
F_periodic = define_forcing_vector(n, Fp_periodic, t)
F_total = F_static + F_periodic

f̂₀ = vec(apply_fft_projection(Eᴴ, F_total))
x₀ = K \ Fp_static
x̂₀ = equilibrium_to_fourier(x₀, n, H)

g, dg = define_nonlinearity(n, [3, 4, 5])

p = (E = E, Eᴴ = Eᴴ, H = H, M = M, C = C, K = K, f̂ = f̂₀, g = g)

cont_pars = ContinuationParameters(λmin = λ₀, λmax = 2.0, Δs = 0.001, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(ϵₓ = 1e-2, ϵᵣ = 1e-6, maxiters = 25), verbose = true)

prob = ContinuationProblem(residual_hbm, cont_pars, p; autodiff = true)
sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ

function amplitud_xi(sol, H, i)
    n = 5
    idx_b1 = n * 1 + i
    idx_a1 = n * 2 + i
    [sqrt(sol.x[idx_b1, j]^2 + sol.x[idx_a1, j]^2) for j in 1:length(sol.λ)]
end

amplitudes_x3 = amplitud_xi(sol, H, 3)

figure = Figure()
ax = Axis(figure[1, 1], limits=(0, 2, 0, 2), xlabel=L"\Omega", ylabel=L"x_3")
lines!(ax, λ_values, amplitudes_x3, color=:dodgerblue, label="Continuation Method")
axislegend(ax, position = :lt)
figure

TFG.save_figure_pdf("scripts/continuation_timeint_resonance.pdf", figure)
