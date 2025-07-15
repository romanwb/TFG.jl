using TFG, ContinuationSuite, GLMakie, DifferentialEquations, ForwardDiff, CairoMakie, LinearAlgebra

" ! NOTAS
Nota 1: Se hace una matriz para guardar los coeficientes de G, de forma que es más comodo partir de esta matriz a la hora de obtener el vector de coeficientes de G ya sea agrupado por nodos o armónicos (solo cambia que el vector se crea cogiendo fila tras fila o columna tras columna)


Nota 2: 

"

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
    #T = eltype(λ) 

    H, Nx, ξ = p.H, p.Nx, p.ξ
    E, Eᴴ = p.E, p.Eᴴ
    γ = p.γ
    # dof para referirse a dof del problema algebraico, no el fisico
    dof_per_node = 2H + 1
    dof_total = Nx * dof_per_node
    Nh = 2H + 1
    #T = eltype(x̂)
    G_mat = zeros(T, Nx, 2H+1) # Nota 1 
    G = zeros(Nx*(2H+1))
    for dof in 1:Nx
        idxs = [dof + (k-1)*Nx for k in 1:Nh]
        coeffs = x̂[idxs]
        G_mat[dof,:] = Eᴴ * (γ*(E * coeffs).^3)
    end
    
    G = vec((G_mat)) #Ordenados por armónicos (coge columna a columna)

    #T = eltype(x̂)
    LHS = zeros(dof_total, dof_total)
    #Kxx diagonal
    LHS[1:Nx, 1:Nx] .= p.Kxx


    for k in 1:H
    #ω² = T.((p.Kaa))
    ω² = p.Kaa
      num1 = ω² .- ((k^2) * λ^2)
      denom = num1.^2 .+ (k * λ * ξ * ω²).^2
      #Dimensiones 2x2
      D1 = Diagonal(num1 ./ denom)
      D2 = Diagonal((k * λ * ξ * ω²) ./ denom)

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

#Params
N, H = 2^9, 10
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
    #T = eltype(λ)
    #Max_prev = p.Max
    #Max = T.(Max_prev)
    Max = p.Max
    H, Nx, ξ = p.H, p.Nx, p.ξ
    γ = p.γ
    Fₓ = zeros(Nx*(2H+1))
    Fₐ = zeros((n-Nx)*(2H+1))

    Fₓ[7] = 3.0
    Fₓ[13] = 3.0

    fₐ = [-1.3763819204711734, 0.3249196962329064]
    Fₐ[[5,6]] .= fₐ
    Fₐ[[9,10]] .= fₐ


    #ω² = T.((p.Kaa))
    ω² = p.Kaa
    F = zeros(Nx*(2H+1))
    for k in 1:H
        num1 = ω² .- (k^2) * (λ^2)
        denom = num1.^2 .+ (k * λ * ω² * ξ).^2
    
        #Dimensiones 2x2
        D1 = Diagonal(num1 ./ denom)
        D2 = Diagonal((k * λ * ω² * ξ) ./ denom)
    
        F[(Nx*2k-(Nx-1)):(Nx*2k)] = Fₓ[(Nx*2k-(Nx-1)):Nx*2k] + (k^2)*(λ^2)*Max'*(D1*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] - D2*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
        F[(Nx*2k+1):(Nx*2k+Nx)] = Fₓ[(Nx*2k+1):(Nx*2k+Nx)] + (k^2)*(λ^2)*Max'*(D2*Fₐ[((n-Nx)*2k-((n-Nx)-1)):((n-Nx)*2k)] - D1*Fₐ[((n-Nx)*2k+1):((n-Nx)*2k+(n-Nx))])
    end    

    return F
end


E, Eᴴ = fft_matrices(N, H)
T = Float64
F= zeros(T, Nx*(2H+1))

p = HBMParams{T}(Kxx, Mxx, Max, Kaa, F, γ, H, Nx, E, Eᴴ, ξ)

x̂₀ = zeros(Nx*(2H+1))

#Solver
λ₀ = 0.25
cont_pars = ContinuationParameters(λmin = λ₀, λmax = 1.75, Δs = 0.01, maxsteps = 10_000,
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
Ωs10new, ampsH10new = extract_amplitude_vs_frequency_by_harmonic_order(sol, E, H, dof, Nx)

begin
fig = Figure()
axis = Axis(fig[1, 1], limits=(0.25, 1.75, 0, 2), xlabel = L"\Omega", ylabel = L"\mathrm{max}(|x_3(t)|)")
#lines!(axis, Ωs3, ampsH3; color = :blue, label = L"H=3")
#lines!(axis, Ωs5, ampsH5; color = :red, label = L"H=5")
#lines!(axis, Ωs10, ampsH10; color = :orange, label = L"H=10")
lines!(axis, Ωs10new, ampsH10new; color = :blue, label = L"H=10new")

lines!(axis, ω_axis, amplitudes_int; color = :black, linestyle = :dash, linewidth = 1, label = L"Time")

axislegend(axis, position = :lt)
fig
end

using CairoMakie
TFG.save_figure_pdf("scripts/real/5dof_nlproblem/5dof_nlproblem.pdf", fig)






function reconstruct_xt(sol, E, H, dof::Int, Nx::Int, idx_sol::Int)
    Nh = 2H + 1
    x̂ = sol.x[:, idx_sol]  # solución para cierto Ω
    coeffs = [x̂[Nx*(k-1) + dof] for k in 1:Nh]
    x_t = E * coeffs
    return x_t
end

idx_sol = 1175# por ejemplo, paso número 150 de la continuación
Ω = sol.λ[idx_sol]
x_t = reconstruct_xt(sol, E, H, 3, Nx, idx_sol)
Nt = length(x_t)
t_vals = LinRange(0, 2π, Nt + 1)[1:end-1]  # para evitar repetir el punto final

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"t", ylabel = L"x_3(t)",
          title = "Respuesta temporal reconstruida (índice = $idx_sol)")
lines!(ax, t_vals, abs.(x_t))
fig
