using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
use_formal_theme!()

@load "scripts/tbm_cyclic/1.first_section/1. problem/H5" Ω Xⱼ Aⱼ p
H = p.H
Nc = 49
op = ProjectorMaker(2^10, p.H)
x_max, a_max= make_FRF(op, Xⱼ, Aⱼ)

a_max[1][1]

# Devuelve las posiciones en A de los valores dados en b: idxsΩ = make_indexs(A, b)
Ωs_to_eval = [0.90, 0.917, 0.92, 0.93, 0.94]
idxsΩ = make_indexs(Ω, Ωs_to_eval)

# Figura con puntos de la curva de interés
begin
    fig = Figure()
    ax = Axis(fig[1,1], limits = (0.89,0.95, nothing,nothing), xlabel = "Ω", ylabel = L"\text{max}|a_{1}|")

    lines!(ax, Ω, a_max[1])
     xs = [Ω[idxsΩ[1]], Ω[idxsΩ[2]], Ω[idxsΩ[3]], Ω[idxsΩ[4]], Ω[idxsΩ[5]]]
     ys = [a_max[1][idxsΩ[1]], a_max[1][idxsΩ[2]], a_max[1][idxsΩ[3]], a_max[1][idxsΩ[4]], a_max[1][idxsΩ[5]]]
     labels = ["$(Ωs_to_eval[1])", "$(Ωs_to_eval[2])", "$(Ωs_to_eval[3])", "$(Ωs_to_eval[4])", "$(Ωs_to_eval[5])",]
     scatter!(ax, xs, ys, color=:white)
     text!(ax, xs, ys; text=labels, align=(:center, :center), offset = (-40, 0))
    fig
end

# using CairoMakie
# save("scripts/tbm_cyclic/1.first_section/3. study_of_contact_plane/curve_with_freq_points.pdf", fig; backend=CairoMakie)



# Devuelve una matriz [filas=Nx, columnas=idxsΩ] con true o false, si está en régimen stick o slip respectivamente
isstick = stick_or_slip(Xⱼ, idxsΩ, Nc, 5, p.kₙ, p.kₜ₁, p.kₜ₂, p.μ, p.xp, H)

# Malla general
f = full_mesh()

# Colorea de negro o rojo los puntos de la malla general

plot_slip_stick(isstick, 1, 7)
using CairoMakie
save("scripts/tbm_cyclic/1.first_section/3. study_of_contact_plane/omega1.pdf", f; backend=CairoMakie)

# Para ver los ciclos de las fuerzas en funcion de x para un nodo Nₓ en concreto
Nₓ = 17
t, x₁, x₂, t_force1, t_force2, w1_vec, w2_vec = study_1contact(Xⱼ, idxsΩ[1], Nₓ, 5, p.kₙ, p.kₜ₁, p.kₜ₂, p.μ, p.xp)

lines(x₁, t_force1)
lines(x₂, t_force2)
lines(t, x₁)
lines(t, t_force1)
lines(t, x₂)
lines(t, t_force2)

lines(t, w1_vec)
w1_vec == zeros(size(w1_vec))

lines(t, w2_vec)
w2_vec == zeros(size(w1_vec))



