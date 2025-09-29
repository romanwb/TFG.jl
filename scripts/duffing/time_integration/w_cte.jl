# Codigo de integración temporal para una w=cte. Evolucion de la amplitud en t
using TFG, DifferentialEquations, GLMakie, CairoMakie

ξ, ϵ, f₀ = 1, 0.01, 1
F₀ = ϵ * f₀
ξ̃ = ϵ * ξ
NPeriods = 30
ω = 1.01

function DuffingOscillator(dy, y, p, t)
    ξ̃, ϵ, F₀, ω = p
    dy[1] = y[2] 
    dy[2] = - ξ̃*y[2] - y[1] - ϵ*y[1]^3 + 2*F₀*cos(ω*t)
end

y0 = [0, 0]
t = NPeriods * 2π / ω
tspan = (0.0, t)
p = (ξ̃, ϵ, F₀, ω)
prob = ODEProblem(DuffingOscillator, y0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.01)
ScalePeriod = NPeriods * sol.t / t

use_formal_theme!()
begin
    f = Figure()
    ax = Axis(f[1, 1], 
    limits = (0, NPeriods, -2, 2),
    xlabel = L"\text{Periodos [rad/s]}", 
    ylabel = L"x(t)", 
    )

    lines!(ax, ScalePeriod, sol[1, :])
    f
end

# -------- Guardar como PDF vectorial --------
CairoMakie.activate!()
save("figures/duffing/time_integration/w_cte/respuesta_duffing_$(NPeriods)_periods.pdf", f; backend=CairoMakie)