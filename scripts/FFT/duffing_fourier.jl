using TFG, ContinuationSuite, GLMakie

N, H = 2^6, 2
Ω, ξ, ϵ = 0.98, 0.01, 0.01
E, Eᴴ = fft_matrices(N, H)

x̂₀ = [0.0, 0.1, 0.1, 0.1, 0.1]
f̂ = [0.0, 0.02, 0.0, 0.0, 0.0]

g(x) = x^3
dg(x) = 3x^2


" Resolucion del problema para una frecuencia Ω dada"

p = DuffingParams(N, H, Ω, ξ, ϵ, E, Eᴴ, x̂₀, f̂, g, dg)
prob = NonlinearProblem(duffing, p)
#prob = NonlinearProblem(duffing, p; jac=jacobian_matrix)
sol_fourier = solve(prob, Newton(), x̂₀)
sol_fisica = E * sol_fourier
 
t = (0:63) .* (2π / 64) #ejeX
f = Figure()
ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")
lines!(ax, t, sol_fisica)
f


" Resolución del problema en un rango de frecuencias: "

resonance_newton=[]
for i in -0.04:0.001:0.04
    p_Ωvar = DuffingParams(N=N, H=H, Ω= 1+i, ξ=ξ, ϵ=ϵ, E=E, Eᴴ=Eᴴ, x₀=x̂₀, f̂=f̂, g=g, dg=dg)
    prob = NonlinearProblem(duffing, p_Ωvar)
    sol = solve(prob, Newton(), x̂₀)
    max_sol_fisica = maximum(E * sol)
    push!(resonance_newton, max_sol_fisica)
end

freqs = (-0.04:0.001:0.04) .+ 1 #ejeX
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")
lines!(ax, freqs, resonance_newton)
fig








# Ejemplo ilustrativo del funcionamiento de la matriz E, INCLUIR ESTO EN TEST
sin_vector = zeros(N)
for i in 1:N
    sin_vector[i] = sin(2π*(i-1)/N)
end

 Eᴴ * sin_vector
