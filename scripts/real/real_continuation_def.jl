using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

#Build T
function build_T(Max, Kaa, Ω, H, ξ)
    T = eltype(Ω)
    Nx = size(Max, 2)
    N_total = (2H + 1) * Nx
    block = zeros(T, N_total, N_total)

    ω² = diag(Kaa)

    for k in 1:H
        num1 = ω² .- T(k^2) * Ω^2
        denom = num1.^2 .+ (T(k) * Ω * ω² * ξ).^2
        #denom .= denom .+ T(1e-12)

        D1 = Diagonal(num1 ./ denom)
        D2 = Diagonal((T(k) * Ω * ω² * ξ) ./ denom)

        Mk1 = -T(k^4) * Ω^4 * Max' * D1 * Max
        Mk2 =  T(k^4) * Ω^4 * Max' * D2 * Max

        T11 = Mk1
        T12 = Mk2
        T21 = -Mk2
        T22 = Mk1

        row_start = (2k - 1) * Nx + 1
        row_end   = (2k + 1) * Nx

        block[row_start:row_end, row_start:row_end] = [T11 T12; T21 T22]
    end

    return block
end



struct HBMParams
    Kxx::Matrix{Float64}
    Mxx::Matrix{Float64}
    Max::Matrix{Float64}
    Kaa::Diagonal{Float64, Vector{Float64}}
    F::Vector{Float64}
    γ::Float64
    H::Int
    Nx::Int
    E::Matrix{Float64}
    Eᴴ::Matrix{Float64}
    ξ::Float64
end


function continuation_system(x̂, λ, p::HBMParams)
    H, Nx, ξ = p.H, p.Nx, p.ξ
    E, Eᴴ = p.E, p.Eᴴ
    γ = p.γ
    dof_per_node = 2H + 1
    dof_total = Nx * dof_per_node

    X_matrix = reshape(x̂, dof_per_node, Nx)'       # (Nx × dof_per_node)
    X_time = X_matrix * E'                          # (Nx × Nt)
    g_time = γ .* X_time.^3
    G = vec(g_time * Eᴴ')  # tamaño Nx*(2H+1)

    
    T = eltype(x̂)
    LHS = zeros(T, dof_total, dof_total)
    #Kxx diagonal
    LHS[1:Nx, 1:Nx] .= p.Kxx


#         ⎡ 1 - k^2*λ^2       kλξ   ⎤
#   A_k = ⎢                         ⎥
#         ⎣     -kλξ     1 - k^2*λ^2⎦
    A = system_matrix(H, ξ, λ)  # 2H × 2H
    I_H = Matrix{T}(I, 2H, 2H)
# Generacion de diagonal de matrices
    LHS_k = kron(A, p.Kxx) + kron(I_H, -λ^2 * p.Mxx)
# NOTA: Recordar que las primeras Nx filas y columnas corresponden al armónico constante k=0
    LHS[Nx+1:end, Nx+1:end] .= LHS_k

    LHS .+= build_T(p.Max, p.Kaa, λ, H, ξ)

    return LHS * x̂ + G - p.F
end

#Params
N, H = 2^6, 3
ξ, ϵ = 0.05, 1
ξ̃ = ϵ * ξ
γ = 1
E, Eᴴ = fft_matrices(N, H)
n = 5
Nh = 2H + 1
m = 1.0
k = 10.0
E, Eᴴ = fft_matrices(2H+1, H)
Nx = 3

# Example_HBM
omega_a2 = [3.8196601125010514, 26.18033988749895]  # → Kaa = Diagonal(omega_a2.^2)
Kaa = Diagonal(omega_a2.^2)

Mxx = [3.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0]

Max = [-1.3763819204711734 0.0 0.0;
        0.3249196962329064 0.0 0.0]

Kxx = [10.0 -10.0 0.0;
       -10.0 20.0 -10.0;
         0.0 -10.0 20.0]

Psi = [1.0 0.0 0.0 1.0 0.0 0.0]

Phi = [-0.850651 -0.525731;
       -0.525731 0.850651]

fa = [-1.3763819204711734, 0.3249196962329064]
fx = [3.0, 0.0, 0.0]
Rx = [-1.2, -0.1, -0.1]
xe0 = zeros(2)

function build_time_domain_force(Nx, Nt)
    t = range(0, 2π, length=Nt+1)[1:end-1]
    fbase = sin.(t) .+ sin.(2t) 

    fx = [3.0, 0.0, 0.0]
    f_time = fx .* fbase'

    return f_time, t
end

function time_to_fourier_force(f_time, Eᴴ)
    # f_time: Nx × Nt
    Fhat_matrix = f_time * Eᴴ'
    Fhat_vector = vec(Fhat_matrix)
    return Fhat_vector
end

function reorder_force_by_harmonic(Fhat_by_node::Vector, Nx::Int, H::Int)
    dof_per_node = 2H + 1
    F_matrix = reshape(Fhat_by_node, Nx, dof_per_node)
    F_reordered = []
    for k in 1:dof_per_node
        push!(F_reordered, F_matrix[:, k]) 
    end

    return vcat(F_reordered...)
end

Nt = 300
E, Eᴴ = fft_matrices(Nt, H)
f_time, t = build_time_domain_force(Nx, Nt)
Fhat_by_node = time_to_fourier_force(f_time, Eᴴ)
F = reorder_force_by_harmonic(Fhat_by_node, Nx, H)
p = HBMParams(Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)

E, Eᴴ = fft_matrices(N, H)

###################################################################
#####################       PRELOAD        ########################
###################################################################

function solve_static_preload(Kxx, Rx, γ; tol=1e-10, maxiter=50)
    x = zeros(length(Rx))  # inicialización
    for iter in 1:maxiter
        gx = γ .* x.^3
        R = Kxx * x + gx + Rx
        if norm(R) < tol
            return x
        end
        # Jacobiano: J = Kxx + diag(3γx²)
        J = Kxx + Diagonal(3γ .* x.^2)
        Δx = -J \ R
        x += Δx
    end
    error("No converge el estado de precarga")
end


function initial_guess_from_preload(Kxx, Rx, γ, H, Nx)
    x_p = solve_static_preload(Kxx, Rx, γ)
    dof_per_node = 2H + 1
    x̂₀ = zeros(Nx * dof_per_node)
    for j in 1:Nx
        x̂₀[(j-1)*dof_per_node + 1] = x_p[j]  # solo el coef. constante
    end
    return x̂₀
end


x̂₀ = initial_guess_from_preload(Kxx, Rx, γ, H, Nx)

#Solver
λ₀ = 0.00
cont_pars = ContinuationParameters(λmin = λ₀, λmax = 2.0, Δs = 0.01, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

prob = ContinuationProblem(continuation_system, cont_pars, p; autodiff = true)

sol = continuation(prob, x̂₀, λ₀)

λ_values = sol.λ
size(sol.λ)
size(sol.x)

function extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof::Int, Nx::Int)
    Nh = 2H + 1
    Nt = size(E, 1)
    Nsteps = length(sol.λ)

    amplitudes = zeros(Nsteps)

    for i in 1:Nsteps
        x̂ = sol.x[:, i]
 
        coeffs = [x̂[Nx*(k-1) + dof] for k in 1:Nh]
        x_t = E * coeffs
        amplitudes[i] = maximum(abs.(x_t))
    end

    freqs = sol.λ
    return freqs, amplitudes
end


function plot_FRF(sol, E, H, dof::Int)
    Ωs, amps = extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof, Nx)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\Omega", ylabel = L"\mathrm{max}(|x(t)|)")
    lines!(ax, Ωs, amps)
    return fig
end


#Plot
fig = plot_FRF(sol, E, H, 1)
display(fig)

using CairoMakie
TFG.save_figure_pdf("scripts/REAL.pdf", fig)
