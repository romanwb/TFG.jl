# Cálculo de la curva de resonancia haciendo uso de integracion temporal. ESTA VERSION ACTUALMENTE NO FUNCIONA, PERO INTEGRA LAS FUNCIONES MEJOR. PENDIENTE DE REVISION


using TFG, DifferentialEquations, GLMakie
#Parametros del problema de Duffing necesarios
ξ, ϵ = 1, 0.01
ξ̃ = ϵ * ξ
H = 1
f̂₀ = zeros(2H+1)
f̂₀[2] = 0.02
g(x) = x^3
λ₀ = 0.96
n = 200 #nº de frecuencias evaluadas en la integracion temporal
y₀ = [0, 0] # CI del problema temporal

function time_integration_values_example(ξ, ϵ, n, λ₀, f̂₀, y₀, g; uniform=true)
    tspan = (0.0, 600.0)
    diff, tol = Inf, 1e-2
    time_integration = zeros(n)
    Δω_axis = zeros(n)

    # Creamos el vector Δω_axis según el tipo de espaciado
    if uniform
        Δω_axis .= LinRange(-4, 4, n)
    else
        for i in 1:n
            Δω_axis[i]=((1-(4e-2)*cos(0.5*π*(2*(i-1)/(n-1))))-1)/ϵ
        end
    end

    for i ∈ 1:n
        old_max = 0
        # Seguimos usando tu constructor con i, n, λ₀, etc.
        p = DuffingParamsTimeIntegration(ξ, ϵ, f̂₀, i, n, λ₀, g)
        tspan = (0.0, 600.0)
        tol = 1e-4
        diff = Inf
        while diff > tol
            step = 25
            prob = ODEProblem(duffing_time_domain, y₀, tspan, p)
            sol = DifferentialEquations.solve(prob, Tsit5(), saveat=0.1)
            new_max = maximum((sol[1, end-step:end]))
            diff, old_max = abs(new_max - old_max), new_max
            tspan = (tspan[1], tspan[2] + step)
        end

        time_integration[i] = old_max
    end

    return Δω_axis, time_integration
end

Δω_axis, time_integration = time_integration_values_example(ξ, ϵ, n, λ₀, f̂₀, y₀, g; uniform=false)
ω_axis  = 1 .+ 0.01 .* Δω_axis

begin
fig = Figure()
ax = Axis(fig[1, 1],   
limits=(0.96, 1.04, 0, 2),
xlabel = L"\text{Frecuencia}", 
ylabel = L"x(t)", 
# title = L"\text{Oscilador de Duffing}",
xticklabelfont = "Latin Modern Roman",
yticklabelfont = "Latin Modern Roman",
)
scatter!(ax, ω_axis, time_integration;
    marker = :circle,
    color = :black,  
    strokecolor = :black,      
    strokewidth = 0.5,           
    markersize = 3,
    label = "Time integration"
)
fig
end

TFG.save_figure_pdf("scripts/time_integration/w_var/respuesta_duffing_w_var.pdf", fig)