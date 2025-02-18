ξ, ϵ, f₀ = 1, 0.01, 2 #try f₀=1 for greater concordance
F₀ = ϵ * f₀
ξ̃ = ϵ * ξ
NPeriods = 500
ω = 1.01
    
#Duffing Oscillator: ẍ + ξ̇x + x + ϵx³ = 2F₀cos(wt) where y[2]=̇x, y[1]=x
function DuffingOscillator(dy, y, p, t)
    ξ̃, ϵ, F₀, ω = p
    dy[1] = y[2] 
    dy[2] = - ξ̃*y[2] - y[1] - ϵ*y[1]^3 + 2*F₀*cos(ω*t)
end
    
y0 = [0, 0] # Initial conditions
t = NPeriods * 2π/ω
tspan= (0.0, t)
p = (ξ̃, ϵ, F₀, ω)
    
prob = ODEProblem(DuffingOscillator, y0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.01)
ScalePeriod = NPeriods * sol.t / t 
    
f = Figure()
ax = Axis(f[1, 1], xlabel = "Periods", ylabel = "𝑥")
lines!(ax, ScalePeriod, sol[1,:])
f