using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

" ! NOTAS
Nota 1: Se hace una matriz para guardar los coeficientes de G, de forma que es más comodo partir de esta matriz a la hora de obtener el vector de coeficientes de G ya sea agrupado por nodos o armónicos (solo cambia que el vector se crea cogiendo fila tras fila o columna tras columna)


Nota 2: 

"


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
    # dof para referirse a dof del problema algebraico, no el fisico
    dof_per_node = 2H + 1
    dof_total = Nx * dof_per_node
    Nh = 2H + 1
    T = eltype(x̂)
    G_mat = zeros(T, Nx, 2H+1) # Nota 1 
    G = zeros(Nx*(2H+1))
    for dof in 1:Nx
        idxs = [dof + (k-1)*Nx for k in 1:Nh]
        coeffs = x̂[idxs]
        G_mat[dof,:] = Eᴴ * ((E * coeffs).^3)
    end
    
    G = vec((G_mat)) #Ordenados por armónicos (coge columna a columna)

    T = eltype(x̂)
    LHS = zeros(T, dof_total, dof_total)
    #Kxx diagonal
    LHS[1:Nx, 1:Nx] .= p.Kxx


    for k in 1:H
    ω² = T.(diag(p.Kaa))

      num1 = ω² .- T(k^2) * λ^2
      denom = num1.^2 .+ (T(k) * λ * ω² * ξ).^2
      #Dimensiones 2x2
      D1 = Diagonal(num1 ./ denom)
      D2 = Diagonal((T(k) * λ * ω² * ξ) ./ denom)

        B11 = p.Kxx - (k^2)*(λ^2)*p.Mxx - T(k^4) * λ^4 * Max' * D1 * Max
        B12 = k*λ*ξ*p.Kxx + T(k^4) * λ^4 * Max' * D2 * Max
        B21 = - k*λ*ξ*p.Kxx - T(k^4) * λ^4 * Max' * D2 * Max
        B22 = p.Kxx - (k^2)*(λ^2)*p.Mxx - T(k^4) * λ^4 * Max' * D1 * Max
        B = [B11 B12; B21 B22]

        LHS[(Nx*2k-(Nx-1)):(Nx*2k+(Nx)), (Nx*2k-(Nx-1)):(Nx*2k+(Nx))] .= B

    end

    Force = external_force(λ, p::HBMParams)
    return LHS * x̂ + G - Force
end

#Params
N, H = 2^6, 5
ξ, ϵ = 0.05, 1
ξ̃ = ϵ * ξ
γ = 1.0
m = 1.0
k = 10.0
E, Eᴴ = fft_matrices(N, H)
n = 5 
Nx = 3 #! nº de masas con friccion no lineal
Nh = 2H + 1


# Example_HBM
omega_a2 = [3.8196601125010514, 26.18033988749895]
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

function external_force(λ, p::HBMParams)
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
        num1 = ω² .- T(k^2) * λ^2
        denom = num1.^2 .+ (T(k) * λ * ω² * ξ).^2
    
        #Dimensiones 2x2
        D1 = Diagonal(num1 ./ denom)
        D2 = Diagonal((T(k) * λ * ω² * ξ) ./ denom)
    
        F[(Nx*2k-(Nx-1)):(Nx*2k)] = Fₓ[(Nx*2k-(Nx-1)):Nx*2k] + (k^2)*(λ^2)*Max'*(D1*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] - D2*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
        F[(Nx*2k+1):(Nx*2k+Nx)] = Fₓ[(Nx*2k+1):(Nx*2k+Nx)] + (k^2)*(λ^2)*Max'*(D2*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] - D1*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
    end    

    return F
end


E, Eᴴ = fft_matrices(N, H)
F= zeros(N)
T = Float64

p = HBMParams{T}(Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)

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

print(λ_values)

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
Ωs5, ampsH5 = extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof, Nx)

begin
fig = Figure()
axis = Axis(fig[1, 1], limits=(0, 2, 0, 2), xlabel = L"\Omega", ylabel = L"\mathrm{max}(|x_3(t)|)")
lines!(axis, Ωs3, ampsH3; color = :blue, label = L"H=3")
lines!(axis, Ωs5, ampsH5; color = :red, label = L"H=5")
lines!(axis, Ωs10, ampsH10; color = :orange, label = L"H=10")



lines!(axis, ω_axis, amplitudes_int; color = :black, linestyle = :dash, linewidth = 1, label = L"Time")

axislegend(axis, position = :lt)
fig
end

using CairoMakie
TFG.save_figure_pdf("scripts/real/5dof_nlproblem/5dof_nlproblem.pdf", fig)
