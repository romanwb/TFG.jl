using TFG, GLMakie

# RECORDATORIO: comentar sobre el bucle j in 1:2 (para hacer la gráfica habría que guardar informacion desde w=0)
function g_completa(x)
    w = 0.0
    t = range(0, 2π, length = length(x) + 1)[1:(end-1)]
    N = ones(length(x))
    N = 1.0 .+ 1.25 .* sin.(t) # N(t) ≠ cte
    # N = 1.0 .+ t
    #  N .*= 2.0 # N(t) = cte
    μ = 0.5
    kₜ = 100.0


    t_full = zeros(2*length(x)) # Recoge historico completo
    t_force = zeros(length(x))

    idx = 1
    for j in 1:2
        for i in eachindex(x)
            t_force[i], w = tangencial_force(x[i], w, kₜ, μ, N[i])
            t_full[idx] = t_force[i]
            idx += 1
        end
    end

    return t_full
end

T = g_completa(x_full)

using JLD2
# @save "scripts/friction_forces/friction_forces.jld2" T_kₜ2 T_kₜ1 T_kₜ5 T_kₜ100
abs_disp = 1
x_1 = range(0, abs_disp, length=100)
x_2 = range(abs_disp, -abs_disp, length=100)
x_3 = range(-abs_disp, 0, length=100)
x_full = vcat(collect(x_1), collect(x_2), collect(x_3))
@load "scripts/friction_forces/friction_forces_N_sin.jld2" T_kₜ2 T_kₜ1 T_kₜ5 T_kₜ100
#  @load "scripts/friction_forces/friction_forces_N_cte.jld2" T_kₜ2 T_kₜ1 T_kₜ5 T_kₜ100

begin
use_formal_theme!()
fig = Figure(resolution = (1200, 1000))
ax1  = Axis(fig[1,1], limits=(-1.2, 1.2, -1.2, 1.4), xlabel=L"x_{t}", ylabel=L"T")
ax2  = Axis(fig[1,2], limits=(-1.2, 1.2, -1.2, 1.4), xlabel=L"x_{t}", ylabel=L"T")
ax3  = Axis(fig[2,1], limits=(-1.2, 1.2, -1.2, 1.4), xlabel=L"x_{t}", ylabel=L"T")
ax4  = Axis(fig[2,2], limits=(-1.2, 1.2, -1.2, 1.4), xlabel=L"x_{t}", ylabel=L"T")


lines!(ax1, x_full, T_kₜ2[1:300]; color=:blue, label=L"\text{k}_{\text{t}} = 2")
lines!(ax1, x_full, T_kₜ2[301:600]; color=:blue)
vlines!(ax1, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)
hlines!(ax1, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)


lines!(ax2, x_full, T_kₜ1[1:300]; color=:green, label=L"\text{k}_{\text{t}} = 1")
lines!(ax2, x_full, T_kₜ1[301:600]; color=:green)
vlines!(ax2, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)
hlines!(ax2, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)

lines!(ax3, x_full, T_kₜ5[1:300]; color=:red, label=L"\text{k}_{\text{t}} = 5")
lines!(ax3, x_full, T_kₜ5[301:600]; color=:red)
vlines!(ax3, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)
hlines!(ax3, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)

lines!(ax4, x_full, T_kₜ100[1:300]; color=:blue, label=L"\text{k}_{\text{t}} = 100")
lines!(ax4, x_full, T_kₜ100[301:600]; color=:blue)
vlines!(ax4, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)
hlines!(ax4, [0.0]; linestyle = :dash, color = :gray, linewidth = 0.8)


axislegend(ax1; position = :lt, labelsize = 24)
axislegend(ax2; position = :lt, labelsize = 24)
axislegend(ax3; position = :lt, labelsize = 24)
axislegend(ax4; position = :lt, labelsize = 24)

display(GLMakie.Screen(), fig)
end

using CairoMakie
save("figures/friction_forces/cycles_n_sin.pdf", fig; backend=CairoMakie)
