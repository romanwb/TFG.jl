using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays, Colors
use_formal_theme!()

using Base: splitpath

function extract_key(path::AbstractString)
    filename = splitpath(path)[end]     # "34x34_ff1.jld2"
    return split(filename, "_")[1]      # "34x34"
end

const Ncs = Dict{String, Int}(
    "7x7"         => 49,
    "34x34"       => 4,
    "345x345"     => 9,
    "23456x23456" => 25,
    "1-7x1-7"     => 4,
    "first_last_cols" => 14,
    "first_last_rows" => 14,
    "first2colums" => 14,
    "upper2rows" => 14
)

function Nc_from_key(key::AbstractString)
    return Ncs[key]
end

path = "scripts/tbm_cyclic/study_of_ext_force/data/EO1/7x7_ff25.jld2"
@load path Ω Xⱼ Aⱼ p
key = extract_key(path)
H, μ, xp  = p.H, p.μ, p.xp
# Nc = 49
Nc = Nc_from_key(key)
op = ProjectorMaker(2^10, p.H)
x_max, a_max= make_FRF(op, Xⱼ, Aⱼ)
lines(Ω, a_max[1])

# Devuelve las posiciones en A de los valores dados en b: idxsΩ = make_indexs(A, b)
# Ωs_to_eval = [0.7, 0.8, 0.85, 0.9, 1.0, 1.1]
Ωs_to_eval = [0.8, 0.83, 0.87, 0.92, 0.96, 1.2]
idxsΩ = make_indexs(Ω, Ωs_to_eval)

# Devuelve una matriz [filas=Nx, columnas=idxsΩ] con true o false, si está en régimen stick o slip respectivamente
n_cycles = 5
# isstick = stick_or_slip(Xⱼ, idxsΩ, Nc, n_cycles, p.kₙ, p.kₜ₁, p.kₜ₂, p.μ, p.xp, H)

t_full, x₁_full, x₂_full, t_force1_full, t_force2_full, w1_vec_full, w2_vec_full, normal_full = study_fullplane(Xⱼ, idxsΩ, Nc, n_cycles, p.kₙ, p.kₜ₁, p.kₜ₂, p.μ, xp, H)

isstick = stick_or_slip2(Nc, idxsΩ, w1_vec_full, w2_vec_full)

normal_ncycles_full = Matrix{Vector{Float64}}(undef, Nc, length(Ωs_to_eval))
for i in 1:Nc
    for j in 1:length(Ωs_to_eval)
    normal_ncycles_full[i,j] = repeat(normal_full[i,j], n_cycles)
    end
end

begin
GLMakie.activate!()
# fig, ax = base_figure(Ω, a_max, idxsΩ, (0.6, 1.2, 50, 700), Ω_labels = string.(Ωs_to_eval))
fig, ax = base_figure(Ω, a_max, idxsΩ, (0.6, 1.2, 0, 500), Ω_labels = string.(Ωs_to_eval))
fig
# Subfiguras
    axs_dots = Vector{Axis}(undef, size(Ωs_to_eval))
    axs_normalF = Vector{Axis}(undef, size(Ωs_to_eval))
    axs_cycle = Vector{Axis}(undef, size(Ωs_to_eval))

    row = [3, 2, 1, 1, 2, 3]
    col = [1, 1, 2, 3, 4, 4]

    sizing = 0.7
    fwidth = sizing*ones(6)
    fheight = sizing*ones(6)
    # fhalign = [1.5, 0.5, 0, 1, 1.5, 1.5]
    # fvalign = [-0.5, 0.5, 0, 0, 0.5, -0.5]
    fhalign = [1,1,0.4,0.4,0.5,0.5]
    fvalign = [0,1,0.2,0.2,1,0]

    k = 1
    for (r, c, w, h, ha, va) in zip(row, col, fwidth, fheight, fhalign, fvalign)
            axs_dots[k] = Axis(fig[r, c], width  = Relative(w), height = Relative(h), halign = ha, valign = va, 
            backgroundcolor = :white, title = "Ω = $(string.(Ωs_to_eval)[k])")
            translate!(axs_dots[k].blockscene, 0, 0, 150)
            full_mesh(axs_dots[k])
            k += 1
    end

    for i in eachindex(Ωs_to_eval)
        # Colorea de negro,rojo o azul los puntos de la malla general

        plot_slip_stick(axs_dots[i], isstick, i, key)
    end

    # Seleccion del contacto a estudiar
    contact_coordinate = 4
    plot_selected_contact.(axs_dots, Ωs_to_eval, contact_coordinate, key)

for i in eachindex(Ωs_to_eval)

    p_start = point_in_fig_from_axis(ax, Ωs_to_eval[i], a_max[1][idxsΩ[i]])
    p_end = center_of_axis_in_fig(axs_dots[i])
    dir = p_end - p_start
    arrows!(fig.scene, [p_start], [dir];
        arrowsize = 20,
        linewidth = 2,
        color = :lightgray,
    )
end

t = t_full[1]
k = 1
    for (r, c, w, h, ha, va) in zip(row, col, fwidth, fheight, fhalign, fvalign)

        if k <= 2
            dha1, dha2 = -1.5, -1.5
            dva1, dva2 = 0.7, -0.7
        elseif k == 3
            dha1, dha2 = -1.3, -1.3
            dva1, dva2 = 0.9, -0.5
        elseif k == 4
            dha1, dha2 = 1.3, 1.3
            dva1, dva2 = 0.9, -0.5
        else # k = 5, 6
            dha1, dha2 = 1.3, 1.3
            dva1, dva2 = 0.7, -0.7
        end

        ha_loc1 = ha + dha1
        ha_loc2 = ha + dha2
        va_loc1 = va + dva1
        va_loc2 = va + dva2


        axs_normalF[k] = Axis(fig[r, c], width  = Relative(0.8*w), height = Relative(0.8*h), halign = ha_loc1, valign = va_loc1, 
        backgroundcolor = :white, xlabel = "t", ylabel = "μN, Fₜ₁")
        lines!(axs_normalF[k], t, μ * normal_ncycles_full[contact_coordinate,k], color=:gray)
        lines!(axs_normalF[k], t, -μ * normal_ncycles_full[contact_coordinate,k], color=:gray)
        lines!(axs_normalF[k], t, t_force1_full[contact_coordinate,k])
        hidedecorations!()

        axs_cycle[k] = Axis(fig[r, c], width  = Relative(0.8*w), height = Relative(0.8*h), halign = ha_loc2, valign = va_loc2, 
        backgroundcolor = :white, xlabel = "t", ylabel = "μN, Fₜ₂")
        lines!(axs_cycle[k], x₁_full[contact_coordinate,k], t_force1_full[contact_coordinate,k])
        hidedecorations!()

            k += 1
    end



display(fig)
end

# using CairoMakie
# save("scripts/tbm_cyclic/1.first_section/3. study_of_contact_plane/full_figure.pdf", fig; backend=CairoMakie)



case_of_study_Ω = 3 # seleccion de uno de los n Ω estudiados
begin
    t = t_full[1]
    # figure = Figure(resolution=(800,800))
    ax = Axis(figure[1,1], xlabel = "t", ylabel = "μN, Fₜ₁")
    lines!(ax, t, μ * normal_ncycles_full[contact_coordinate,case_of_study_Ω], color=:gray)
    lines!(ax, t, -μ * normal_ncycles_full[contact_coordinate,case_of_study_Ω], color=:gray)
    lines!(ax, t, t_force1_full[contact_coordinate,case_of_study_Ω])
    # figure
end

begin
    t = t_full[1]
    figure = Figure(resolution=(800,800))
    ax = Axis(figure[1,1], xlabel = "t", ylabel = "μN, Fₜ₂")
    lines!(ax, t, μ * normal_ncycles_full[contact_coordinate,case_of_study_Ω], color=:gray)
    lines!(ax, t, -μ * normal_ncycles_full[contact_coordinate,case_of_study_Ω], color=:gray)
    lines!(ax, t, t_force2_full[contact_coordinate,case_of_study_Ω])
    figure
end

# begin
    Nc_plot = [39]
    # Nc_plot = [1, 2, 3, 4, 5, 6, 7]

    case_of_study_Ω = 6
    figure2 = Figure()
    ax22 = Axis(figure2[1,1], xlabel = "x₁", ylabel = "t₁")
    for i in Nc_plot
    lines!(ax22, x₁_full[i,case_of_study_Ω], t_force1_full[i,case_of_study_Ω], linewidth=1.8, label = "Nc = $(i)")
    end
    axislegend()
    figure2
# end

# lines(t, w1_vec)
# w1_vec == zeros(size(w1_vec))

# lines(t, w2_vec)
# w2_vec == zeros(size(w1_vec))

# lines(x₂, t_force2)
# # lines(t, x₁)
# lines(t, t_force1)
# # lines(t, x₂)
# lines(t, t_force2)


isstick[:,3]

lines(t, w1_vec_full[3,6])
w1_vec_full[4,6][end-2000:end]

x = (3,4)
eachindex(x)
for i in (x)
    println(i)
end