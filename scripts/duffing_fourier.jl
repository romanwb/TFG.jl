"Solving Duffing Oscillator using Fourier:
    1º Se desarrolla en serie la solucion x(t) (En una serie ∈ Re)
    2º Se sustituye y proyectando sobre las diferentes direcciones (para los diferentes indices k) se monta un sistema de ecuaciones
    3º Haciendo cero los coeficientes del desarrollo de la funcion g(x) no lineal, se consigue una solucion inicial aproximada para los x̂ₖ
    4º Se hace la inversa de estos x̂ₖ para llegar a una aproximada x(t) que se evalua en g(x(t)) obteniendo ésta aproximadamente
    5º Con g(x) → ĝₖ. Ahora se puede iterar este proceso y mejorar los x̂ₖ

    ! Una nota interesante es que, se puede montar facilmente un sistema de k ecuaciones aprovechando la periocidad de las derivadas de funciones periodicas,
      y es esto lo que se hace

"
ξ = 0.01
F₀ = 0.02
Ω = 1.01
ϵ = 0.01
k_max = 5 

# Definidas global para poder usarlas en el prob de ContinuationSuite
global ĝₖ = zeros(k_max + 1)
global ĥₖ = zeros(k_max + 1)


function f(x, p)
    ξ, Ω, F₀ = p

    # x = [a₀, a₁..aₖₘₐₓ, b₀, b₁..bₖₘₐₓ]
    â = x[1 : k_max+1]
    b̂ = x[k_max+2 : 2*(k_max+1)]

    # ! Importante: eqs debe tener el mismo tipo que x[1], 
    # por si ForwardDiff pasa Dual en vez de Float64.
    T = eltype(x)  
    eqs = Vector{T}(undef, 2*(k_max + 1))

    # k=0 => eqs[1] y eqs[k_max+2], típicamente a₀ + g₀=0 y b₀ + h₀=0
    eqs[1] = â[1] + ĝₖ[1]
    eqs[k_max + 2] = b̂[1] + ĥₖ[1]

    # k=1..k_max
    for k in 1:k_max

        # cos(kΩt)
        eqs[k + 1] = - (k^2 * Ω^2)*â[k+1] + (ξ * k * Ω)*b̂[k+1] + â[k+1] + ĝₖ[k+1]
        if k == 1
            eqs[k + 1] -= 2 * F₀
        end

        # sin(kΩt)
        eqs[k_max + k + 2 ] = - (k^2 * Ω^2)*b̂[k+1] - (ξ * k * Ω)*â[k+1] + b̂[k+1] + ĥₖ[k+1]
    end

    return eqs
end


p₀ = (ξ, Ω, F₀)
x₀ = zeros(2*(k_max + 1))
prob = ContinuationSuite.NonlinearProblem(f, p₀)

#Solucion con ĝₖ = 0
sol = ContinuationSuite.solve(prob, Newton(), x₀)

#Valores para fft
N       = 1024
T0      = 1 / Ω
T_total = 200 * T0
time    = range(0, T_total, length=N)

#Las funciones de FFTW no trabajan con reales y es enredado usar éstas, mas simple:
function reconstruct_x_t(x_vec, t, k_max)
    â = x_vec[1 : k_max+1]
    b̂ = x_vec[k_max+2 : 2*(k_max+1)]
    x_t = zeros(length(t))
    for i in eachindex(t)
        τ = t[i]
        tmp = zero(eltype(x_vec))
        for k in 0:k_max
            tmp += â[k+1] * cos(k*Ω*τ) + b̂[k+1] * sin(k*Ω*τ)
        end
        x_t[i] = tmp
    end
    return x_t
end

function extract_coeffs(g_t, k_max)
    G = fft(g_t) ./ length(g_t)
    ghat = real.(G[1 : k_max+1])
    hhat = imag.(G[1 : k_max+1])
    return ghat, hhat
end

# Iteraciones del problema para conseguir una convergencia
tol         = 1e-6
max_iters   = 50
min_iters   = 10
iter        = 0
convergencia = false

while iter < max_iters && (!convergencia || iter < min_iters)
    iter += 1
    x_t = reconstruct_x_t(sol, time, k_max)

    # Def de la funcion no lineal
    g_t = ϵ .* (x_t.^3)

    global ĝₖ, ĥₖ
    ĝₖ, ĥₖ = extract_coeffs(g_t, k_max)
    sol_nueva = ContinuationSuite.solve(prob, Newton(), sol)

    error = norm(sol_nueva .- sol)
    println("Iter $iter: Error = $error")
    if error < tol && iter >= min_iters
        convergencia = true
        println("Convergencia en iter $iter (error = $error).")
    end

    sol = sol_nueva
end


x_t_final = reconstruct_x_t(sol, time, k_max)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time [s]", ylabel="x(t)", title="Duffing (Balance Armónico, con ForwardDiff)")
lines!(ax, time, x_t_final)
fig
