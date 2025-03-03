using DifferentialEquations, CairoMakie

ξ, ϵ, f₀ = 1, 0.01, 1 #try f₀=1 for greater concordance
F₀ = ϵ * f₀
ξ̃ = ϵ * ξ
    
#Duffing Oscillator: ẍ + ξ̇x + x + ϵx³ = 2F₀cos(wt) where y[2]=̇x, y[1]=x
function DuffingOscillator(dy, y, p, t)
    ξ̃, ϵ, F₀, ω = p
    dy[1] = y[2] 
    dy[2] = - ξ̃*y[2] - y[1] - ϵ*y[1]^3 + 2*F₀*cos(ω*t)
end

y0 = [0, 0] # C.I.
tspan = (0.0, 300.0)
diff, tol = Inf, 1e-2
n = 41 #nº de frecuencias evaluadas
old_max_values = zeros(n)
for i ∈ 1:n
    old_max = 0
    p = ( ξ̃, ϵ, F₀, 1-(4e-2)*cos(0.5*π*(2*(i-1)/(n-1))) )
    tspan = (0.0, 300.0)
    tol = 1e-4
    diff = Inf
    while diff > tol
        step = 25
        prob = ODEProblem(DuffingOscillator, y0, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.1)
        new_max = maximum(sol[1, end-step:end])
        diff, old_max = abs(new_max - old_max), new_max
        tspan = (tspan[1], tspan[2] + step)
    end
    old_max_values[i]=old_max
end



# Curva de la solución asintótica ()
Avals = range(1e-3, 2.2, length=3000)
deltaOmega_plus  = Float64[]
deltaOmega_minus = Float64[]

for A in Avals
    radicando = (f₀/(2*A))^2 - (ξ/2)^2
    if radicando >= 0
        raiz = sqrt(radicando)
        push!(deltaOmega_plus,  (3/2)*A^2 + raiz)
        push!(deltaOmega_minus, (3/2)*A^2 - raiz)
    else
        push!(deltaOmega_plus,  NaN)
        push!(deltaOmega_minus, NaN)
    end
end

# Plot
# Se divide entre ϵ, ya que ωf = 1 + ϵ∆ωf 
# Aqui se supone que la figura A.2 debería tener un epsilon en el eje x: ϵ∆ωf
Δω_axis=zeros(n)
    for i in 1:n
        Δω_axis[i]=((1-(4e-2)*cos(0.5*π*(2*(i-1)/(n-1))))-1)/ϵ
    end
f = Figure()
ax = Axis(f[1, 1], limits=(-4, 4, nothing, nothing), xlabel = "Δω", ylabel = "2|A|")
lines!(ax, deltaOmega_plus, [2*A for A in Avals], label="Rama +")
lines!(ax, deltaOmega_minus, [2*A for A in Avals], label="Rama -")
scatter!(Δω_axis, old_max_values)
f

println(Δω_axis[25])
println(old_max_values[25])