# Recibe el idx del contacto (idxNᵪ) y los idx del vector Ω correspondiente a las frecuencias a estudiar (idxΩ). También otros valores del problema
function study_1contact(Xⱼ, idxΩ, idxNᵪ, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
    X₁ = Xⱼ[(1 + (2H + 1) * 3(idxNᵪ - 1)):(1 + 2H + (2H + 1) * 3(idxNᵪ - 1)), :]
    X₂ = Xⱼ[(2H + 2 + (2H + 1) * 3(idxNᵪ - 1)):(4H + 2 + (2H + 1) * 3(idxNᵪ - 1)), :]
    X₃ = Xⱼ[(4H + 3 + (2H + 1) * 3(idxNᵪ - 1)):(6H + 3 + (2H + 1) * 3(idxNᵪ - 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, 5)
    # xp = hbm_prob.p.preload_state.xₚ
    x₁ = (E * X₁)[:, idxΩ] .+ xp[1 + 3 * (idxNᵪ - 1)]
    x₂ = (E * X₂)[:, idxΩ] .+ xp[2 + 3 * (idxNᵪ - 1)]
    x₃ = (E * X₃)[:, idxΩ] .+ xp[3 + 3 * (idxNᵪ - 1)]
    # x₁ = (E * X₁)[:, idxΩ]
    # x₂ = (E * X₂)[:, idxΩ]
    # x₃ = (E * X₃)[:, idxΩ]

    t_force1 = zeros(2^10 * n_cycles)
    t_force2 = zeros(2^10 * n_cycles)
    normal = HBMContinuation.normal_force.(x₃, kₙ, 0.0)
    w1_vec = zeros(2^10 * n_cycles)
    w2_vec = zeros(2^10 * n_cycles)
    ws1, ws2 = 0.0, 0.0
    for j in 1:n_cycles
        for i in eachindex(x₁)
            t_force2[i + 2^10 * (j - 1)], w2 = HBMContinuation.tangential_force(
                x₂[i], ws2, kₜ₂, μ, normal[i])
            ws2 = w2
            w2_vec[i + 2^10 * (j - 1)] = w2

            t_force1[i + 2^10 * (j - 1)], w1 = HBMContinuation.tangential_force(
                x₁[i], ws1, kₜ₁, μ, normal[i])
            ws1 = w1
            w1_vec[i + 2^10 * (j - 1)] = w1
        end
    end
    t = range(0, stop = 2π * n_cycles, length = 2^10 * n_cycles)

    x₁ = repeat(x₁, n_cycles)
    x₂ = repeat(x₂, n_cycles)
    return t, x₁, x₂, t_force1, t_force2, w1_vec, w2_vec
end

# Estudia al completo el plano de contacto, para todos los Nx 
function study_fullplane(Xⱼ, idxsΩ, Nc, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
    m, n = Nc, length(idxsΩ)
    t_full = Matrix{Vector{Float64}}(undef, m, n)
    x₁_full = Matrix{Vector{Float64}}(undef, m, n)
    x₂_full = Matrix{Vector{Float64}}(undef, m, n)
    t_force1_full = Matrix{Vector{Float64}}(undef, m, n)
    t_force2_full = Matrix{Vector{Float64}}(undef, m, n)
    w1_vec_full = Matrix{Vector{Float64}}(undef, m, n)
    w2_vec_full = Matrix{Vector{Float64}}(undef, m, n)

    for i in 1:Nc
        for j in eachindex(idxsΩ)
            t_full[i,j], x₁_full[i,j], x₂_full[i,j], t_force1_full[i,j], t_force2_full[i,j], w1_vec_full[i,j], w2_vec_full[i,j] = study_1contact(Xⱼ, idxsΩ[j], i, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
        end
    end
    return t_full, x₁_full, x₂_full, t_force1_full, t_force2_full, w1_vec_full, w2_vec_full
end

function nearest_index(vector, value)
  argmin(abs.(vector .- value))
end

# Recibe los valores de frecuencia que se quieren estudiar y el vector de frecuencias Ω. Da las posiciones en Ω de las frecuencias deseadas.
function make_indexs(vector, vector_of_values)
    idxs = Vector{Int}(undef, length(vector_of_values)) 
    for i in eachindex(vector_of_values)
        idxs[i] = (nearest_index(vector, vector_of_values[i]))
    end
    return idxs
end

function full_mesh()
    f = Figure(;
    figure_padding = 10
    )
    ax = Axis(f[1,1])
    # hidedecorations!(ax)
    x = 1:7
    y = 1:7
    points = Point2f.((x[i], y[j]) for i in eachindex(x) for j in eachindex(y))
    scatter!(ax, points, strokecolor=:white, color=:lightgray, markersize = 15)
    hidedecorations!()
    return f
end

function stick_or_slip(Xⱼ, idxsΩ, Nc, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
    t_full, x₁_full, x₂_full, t_force1_full, t_force2_full, w1_vec_full, w2_vec_full = study_fullplane(Xⱼ, idxsΩ, Nc, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)

    m, n = Nc, length(idxsΩ)
    isstick = Matrix{Bool}(undef, m, n)
    for i in 1:Nc
        for j in eachindex(idxsΩ)
            w1 = w1_vec_full[i,j]
            w2 = w2_vec_full[i,j]
            isstick[i,j] = (all(w1 .== 0.0) && all(w2 .== 0.0))
        end
    end
    return isstick
end

# dim refiere a la dimension del cuadrado, que si son Nc=49 es dim=7
function plot_slip_stick(isstick, n, dim)
    isstick_n = isstick[:,n]
    x = 1:7
    y = 1:7
    for i in eachindex(x)
        for j in eachindex(y)
            if isstick_n[dim*(i-1)+j] == true
                scatter!((i,j), strokecolor=:white, color=:black, markersize = 15)
            else
                scatter!((i,j), strokecolor=:white, color=:red, markersize = 15)
            end
        end
    end
end