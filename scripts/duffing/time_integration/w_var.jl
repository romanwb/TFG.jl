# Cálculo de la curva de resonancia haciendo uso de integracion temporal.
using TFG, DifferentialEquations, GLMakie, JLD2, FileIO, Distributions

#Parametros del problema de Duffing necesarios
p = (
    ξ = 1.0,
    ϵ = 0.01,
    f̂₀ = [0.0, 0.02, 0.0, 0.0, 0.0],
    λ₀ = 0.96,
    n = 60, #nº de frecuencias evaluadas en la integracion temporal
    y₀ = [0, 0] # CI del problema temporal
    )
g(x) = x^3


# Definicion de la distribucion de puntos: tipo normal
function grid_truncnormal(a, b, N; μ = (a+b)/2, σ = (b-a)/6, incluir_extremos=false)
    d = truncated(Normal(μ, σ), a, b)
    p = incluir_extremos ? range(0, 1; length=N) : ((1:N) .- 0.5) ./ N
    return quantile.(Ref(d), p)
end
Δω_axis = grid_truncnormal(p.λ₀, 1+(1-p.λ₀), p.n; σ=(10-0)/8, incluir_extremos=true)

time_integration = time_integration_values(p.ξ, p.ϵ, p.n, p.λ₀, p.f̂₀, p.y₀, g, Δω_axis)
ω_axis  = Δω_axis .- 1


#  @save "data/time_integration/w_var_data.jld2" ω_axis time_integration p
@load "data/time_integration/w_var_data.jld2" ω_axis time_integration p

use_formal_theme!()
begin
    fig = Figure()
    ax = Axis(fig[1, 1],   
    limits=(-0.04, 0.04, 0, 2),
    xlabel = L"Δω_{f}", 
    ylabel = L"|x(t)|", 
    # title = L"\text{Oscilador de Duffing}",
    xticklabelfont = "Latin Modern Roman",
    yticklabelfont = "Latin Modern Roman",
    # xtickformat = xs -> string.(round.(xs .- 1; digits = 3)),
    )
    scatter!(ax, ω_axis, time_integration;
        # marker = :circle,
        # color = :black,  
        # strokecolor = :black,      
        # strokewidth = 0.5,           
        #  markersize = 3,
        # label = "Time integration"
         markersize = 7,

    )
    fig
end
using CairoMakie
save("figures/duffing/time_integration/w_var/respuesta_duffing_w_var.pdf", fig; backend=CairoMakie)

