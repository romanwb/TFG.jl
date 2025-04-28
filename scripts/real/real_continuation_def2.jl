using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

#Build T
function build_T(Max, Kaa, Ω, H, ξ)
    T = eltype(Ω)
    Nx = size(Max, 2)
    N_total = (2H + 1) * Nx
    block = zeros(T, N_total, N_total)

    #ω² = diag(Kaa)
    ω² = Kaa
    I_A = Matrix(I, size(Kaa))
    for k in 1:H
        num1 = ω² .- (T(k^2) * Ω^2 * I_A)
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

#struct HBMParams
#    Kxx::Matrix{Float64}
#    Mxx::Matrix{Float64}
#    Max::Nx*(k-1) + dof] for k in 1:Nh]Matrix{Float64}
#    Kaa::Diagonal{Float64, Vector{Float64}}
#    F::Vector{Float64}
#    γ::Float64
#    H::Int
#    Nx::Int
#    E::Matrix{Float64}
#    Eᴴ::Matrix{Float64}
#    ξ::Float64
#end

struct HBMParams{T}
    Kxx::Matrix{T}
    Mxx::Matrix{T}
    Max::Matrix{T}
    Kaa::Diagonal{T, Vector{T}}
    F::Vector{T}
    γ::T
    H::Int
    Nx::Int
    E::Matrix{T}
    Eᴴ::Matrix{T}
    ξ::T
end



function continuation_system(x̂, λ, p::HBMParams)
    T = eltype(λ) 

    H, Nx, ξ = p.H, p.Nx, p.ξ
    E, Eᴴ = p.E, p.Eᴴ
    γ = p.γ
    dof_per_node = 2H + 1
    dof_total = Nx * dof_per_node
    Nh = 2H + 1
    T = eltype(x̂)
    G_mat = zeros(T, Nx, 2H+1)    
    G = zeros(Nx*(2H+1))
    for dof in 1:Nx
        idxs = [dof + (k-1)*Nx for k in 1:Nh]   # ESTA es la forma correcta
        coeffs = x̂[idxs]
        G_mat[dof,:] = Eᴴ * ((E * coeffs).^3)
    end
    
    G = vec((G_mat))
  # G1 = Eᴴ * ((E * x̂[[1, 4, 7, 10, 13, 16, 19]]).^3)
  # G2 = Eᴴ * ((E * x̂[[2, 5, 8, 11, 14, 17, 20]]).^3)
  # G3 = Eᴴ * ((E * x̂[[3, 6, 9, 12, 15, 18, 21]]).^3)
  # G = [G1[1], G2[1], G3[1], G1[2], G2[2], G3[2], G1[3], G2[3], G3[3], 
  #     G1[4], G2[4], G3[4], G1[5], G2[5], G3[5], G1[6], G2[6], G3[6], 
  #     G1[7], G2[7], G3[7]]

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

    
    Force = final_force(λ, p::HBMParams)
    #Fₓ = zeros(T, Nx*(2H+1))
    #Fₓ[7] = 3
    #Fₓ[13] = 3
    return LHS * x̂ + G - Force
end













#Params
N, H = 2^6, 3
ξ, ϵ = 0.05, 1
ξ̃ = ϵ * ξ
γ = 1.0
E, Eᴴ = fft_matrices(N, H)
n = 5
Nh = 2H + 1
m = 1.0
k = 10.0
Nx = 3
size(E)

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

function final_force(λ, p::HBMParams)
    T = eltype(λ)
    Max_prev = p.Max
    Max = T.(Max_prev)
    H, Nx, ξ = p.H, p.Nx, p.ξ
    γ = p.γ
    Fₓ = zeros(T, Nx*(2H+1))
    Fₐ = zeros(T, (n-Nx)*(2H+1))

    Fₓ[7] = 3
    Fₓ[13] = 3

    fₐ = T.([-1.3763819204711734, 0.3249196962329064])
    Fₐ[[5,6]] .= fₐ
    Fₐ[[9,10]] .= fₐ


    ω² = T.(diag(p.Kaa))

    F = zeros(T, Nx*(2H+1))
    for k in 1:H
        num1 = ω² .- T(k^2) * λ₀^2
        denom = num1.^2 .+ (T(k) * λ₀ * ω² * ξ).^2
    
        #Dimensiones 2x2
        D1 = Diagonal(num1 ./ denom)
        D2 = Diagonal((T(k) * λ₀ * ω² * ξ) ./ denom)
    
        F[[Nx*k+1, Nx*k+2, Nx*k+3]] = Fₓ[[Nx*k, Nx*k+1, Nx*k+2]] + (k^2)*(λ^2)*Max'*(D1*Fₐ[[(n-Nx)*k+1, (n-Nx)*k+2]] - D2*Fₐ[[(n-Nx)*k+3, (n-Nx)*k+4]])
        F[[Nx*k+4, Nx*k+5, Nx*k+6]] = Fₓ[[Nx*k+3, Nx*k+4, Nx*k+5]] + (k^2)*(λ^2)*Max'*(D2*Fₐ[[(n-Nx)*k+1, (n-Nx)*k+2]] - D1*Fₐ[[(n-Nx)*k+3, (n-Nx)*k+4]])
    end    

    return F
end


E, Eᴴ = fft_matrices(N, H)
F= zeros(N)
#p = HBMParams(Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)
T = Float64
#p = HBMParams(T.(Kxx), T.(Mxx), T.(Max), Diagonal(T.(diag(Kaa))), zeros(T, Nx*(2H+1)), T(γ), H, Nx, T.(E), T.(Eᴴ), T(ξ))

p = HBMParams{T}(Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)


# Preload
#function solve_static_preload(Kxx, Rx, γ; tol=1e-10, maxiter=50)
#    x = zeros(length(Rx))  # inicialización
#    for iter in 1:maxiter
#        gx = γ .* x.^3bloque
#    error("No converge el estado de precarga")
#end
#
#
#function initial_guess_from_preload(Kxx, Rx, γ, H, Nx)
#    x_p = solve_static_preload(Kxx, Rx, γ)
#    dof_per_node = 2H + 1
#    x̂₀ = zeros(Nx * dof_per_node)
#    for j in 1:Nx
#        x̂₀[(j-1)*dof_per_node + 1] = x_p[j]  # solo el coef. constante
#    end
#    return x̂₀
#end
#
#
#x̂₀ = initial_guess_from_preload(Kxx, Rx, γ, H, Nx)
#for i in eachindex(x̂₀)
#    println(x̂₀[i])
#end

x̂₀ = zeros(Nx*(2H+1))


#Solver
λ₀ = 0.0
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


dof=3
Ωs, amps = extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof, Nx)
fig = Figure()
axis = Axis(fig[1, 1], limits=(0, 2, 0, 2), xlabel = L"\Omega", ylabel = L"\mathrm{max}(|x(t)|)")
lines!(axis, Ωs, amps)
fig

using CairoMakie
TFG.save_figure_pdf("scripts/REAL.pdf", fig)
