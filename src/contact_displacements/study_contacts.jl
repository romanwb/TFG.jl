
# Recibe el idx del contacto (idxNᵪ) y los idx del vector Ω correspondiente a las frecuencias a estudiar (idxΩ). También otros valores del problema
function study_1contact(Xⱼ, idxΩ, idxNᵪ, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
    X₁ = Xⱼ[(1 + (2H + 1) * 3 * (idxNᵪ - 1)):(1 + 2H + (2H + 1) * 3 * (idxNᵪ - 1)), :]
    X₂ = Xⱼ[(2H + 2 + (2H + 1) * 3 * (idxNᵪ - 1)):(4H + 2 + (2H + 1) * 3 * (idxNᵪ - 1)), :]
    X₃ = Xⱼ[(4H + 3 + (2H + 1) * 3 * (idxNᵪ - 1)):(6H + 3 + (2H + 1) * 3 * (idxNᵪ - 1)), :]
    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
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

    x₁_final = repeat(x₁, n_cycles)
    x₂_final = repeat(x₂, n_cycles)
    return t, x₁_final, x₂_final, t_force1, t_force2, w1_vec, w2_vec, normal
end

# Estudia al completo el plano de contacto, para todos los Nx 
function study_fullplane(Xⱼ, idxsΩ, Nc, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
    m, n = Nc, length(idxsΩ)
    t_full = Matrix{Vector{Float64}}(undef, m, n)
    x₁_full = Matrix{Vector{Float64}}(undef, m, n)
    x₂_full = Matrix{Vector{Float64}}(undef, m, n)
    t_force1_full = Matrix{Vector{Float64}}(undef, m, n)
    t_force2_full = Matrix{Vector{Float64}}(undef, m, n)
    normal_full = Matrix{Vector{Float64}}(undef, m, n)
    w1_vec_full = Matrix{Vector{Float64}}(undef, m, n)
    w2_vec_full = Matrix{Vector{Float64}}(undef, m, n)

    for i in 1:Nc
        for j in eachindex(idxsΩ)
            t_full[i,j], x₁_full[i,j], x₂_full[i,j], t_force1_full[i,j], t_force2_full[i,j], w1_vec_full[i,j], w2_vec_full[i,j], normal_full[i,j] = study_1contact(Xⱼ, idxsΩ[j], i, n_cycles, kₙ, kₜ₁, kₜ₂, μ, xp, H)
        end
    end
    return t_full, x₁_full, x₂_full, t_force1_full, t_force2_full, w1_vec_full, w2_vec_full, normal_full
end