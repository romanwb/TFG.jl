using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

" ! NOTAS
Nota 1: Se hace una matriz para guardar los coeficientes de G, de forma que es más comodo partir de esta matriz a la hora de obtener el vector de coeficientes de G ya sea agrupado por nodos o armónicos (solo cambia que el vector se crea cogiendo fila tras fila o columna tras columna)


Nota 2: 

"

struct HBMParams{T}
    Ttype::Type{T}
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
        G_mat[dof,:] = Eᴴ * (γ*(E * coeffs).^3)
    end
    
    G = vec((G_mat)) #Ordenados por armónicos (coge columna a columna)

    T = eltype(x̂)
    LHS = zeros(T, dof_total, dof_total)
    #Kxx diagonal
    LHS[1:Nx, 1:Nx] .= p.Kxx
    I_A = Matrix(I, size(Kaa))

    for k in 1:H
    #ω² = T.((p.Kaa))
    ω² = p.Kaa
      num1 = ω² .- (T(k^2) * λ^2 * I_A)
      denom = num1.^2 .+ ((T(k) * λ * ω² * ξ).^2)
      #Dimensiones 2x2
      D1 = Diagonal(num1 ./ denom)
      D2 = Diagonal((T(k) * λ * ω² * ξ) ./ denom)

        B11 = p.Kxx - (k^2)*(λ^2)*p.Mxx - T(k^4) * (λ^4) * Max' * D1 * Max
        B12 = k*λ*ξ*p.Kxx + T(k^4) * (λ^4) * Max' * D2 * Max
        B21 = - k*λ*ξ*p.Kxx - T(k^4) * (λ^4) * Max' * D2 * Max
        B22 = p.Kxx - (k^2)*(λ^2)*p.Mxx - T(k^4) * (λ^4) * Max' * D1 * Max
        B = [B11 B12; B21 B22]

        LHS[(Nx*2k-(Nx-1)):(Nx*2k+(Nx)), (Nx*2k-(Nx-1)):(Nx*2k+(Nx))] .= B

    end

    Force = external_force(λ, p::HBMParams)
    return LHS * x̂ + G - Force
end

using Printf
function imprimir_matriz(A)
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            @printf("%8.2f ", A[i, j])
        end
        println()  # Salto de línea al final de cada fila
    end
end

#Params
N, H = 2^9, 3
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
#Kaa = Diagonal(omega_a2.^2)
Kaa = Diagonal(omega_a2)

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

    Fₓ[7] = 3.0
    Fₓ[13] = 3.0

    fₐ = T.([-1.3763819204711734, 0.3249196962329064])
    Fₐ[[5,6]] .= fₐ
    Fₐ[[9,10]] .= fₐ


    #ω² = T.((p.Kaa))
    ω² = p.Kaa
    I_A = Matrix(I, size(Kaa))
    F = zeros(T, Nx*(2H+1))
    for k in 1:H
        num1 = ω² .- (T(k^2) * λ^2 * I_A)
        denom = num1.^2 .+ (T(k) * λ * ω² * ξ).^2
    
        #Dimensiones 2x2
        D1 = Diagonal(num1 ./ denom)
        D2 = Diagonal((T(k) * λ * ω² * ξ) ./ denom)
    
        F[(Nx*2k-(Nx-1)):(Nx*2k)] = Fₓ[(Nx*2k-(Nx-1)):Nx*2k] + (k^2)*(λ^2)*Max'*(D1*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] - D2*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
        F[(Nx*2k+1):(Nx*2k+Nx)] = Fₓ[(Nx*2k+1):(Nx*2k+Nx)] + (k^2)*(λ^2)*Max'*(D2*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] + D1*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
    end    

    return F
end


E, Eᴴ = fft_matrices(N, H)
T = Float64
F= zeros(T, Nx*(2H+1))

p = HBMParams{T}(T, Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)

x̂₀ = zeros(Nx*(2H+1))

#Solver
λ₀ = 0.25
cont_pars = ContinuationParameters(λmin = λ₀, λmax = 1.75, Δs = 0.01, maxsteps = 10_000,
    direction = :forward, predictor = PseudoArcLength(),
    corrector = Newton(), verbose = true)

prob = ContinuationProblem(continuation_system, cont_pars, p; autodiff = false)

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
Ωs3, ampsH3 = extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof, Nx)

begin
fig = Figure()
axis = Axis(fig[1, 1], limits=(0.25, 1.75, 0, 2), xlabel = L"\Omega", ylabel = L"\mathrm{max}(|x_3(t)|)")
lines!(axis, Ωs3, ampsH3; color = :blue, label = L"H=3")
lines!(axis, Ωs5, ampsH5; color = :red, label = L"H=5")
lines!(axis, Ωs10, ampsH10; color = :orange, label = L"H=10")
#lines!(axis, Ωs10new, ampsH10new; color = :blue, label = L"H=10new")

lines!(axis, ω_axis, amplitudes_int; color = :black, linestyle = :dash, linewidth = 1, label = L"Time")

axislegend(axis, position = :lt)
fig
end

using CairoMakie
TFG.save_figure_pdf("scripts/real/5dof_nlproblem/5dof_nlproblemv2.pdf", fig)
~n 

function reconstruct_xt_to_match_integration(sol, E, H, dof::Int, Nx::Int, idx_sol::Int, t_integration_max::Float64)
    Nh = 2H + 1
    x̂ = sol.x[:, idx_sol]
    coeffs = [x̂[Nx*(k-1) + dof] for k in 1:Nh]
    x_period = E * coeffs  # 1 periodo reconstruido

    Ω = sol.λ[idx_sol]
    T_HBM = 2π / Ω
    Nt = length(x_period)

    N_cycles = ceil(Int, t_integration_max / T_HBM)
    x_extended = repeat(x_period, N_cycles)
    total_time = N_cycles * T_HBM
    t_vals = LinRange(0, total_time, N_cycles * Nt + 1)[1:end-1]

    return t_vals, x_extended
end

T_max = maximum(solx5)
idx_sol = 1175# por ejemplo, paso número 150 de la continuación
Ω = sol.λ[idx_sol]
t_hbm, x_hbm = reconstruct_xt_to_match_integration(sol, E, H, 3, Nx, idx_sol, T_max)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"t", ylabel = L"x_3(t)",
          title = "Respuesta temporal reconstruida (índice = $idx_sol)")
lines!(ax, t_vals, abs.(x_t))
fig
