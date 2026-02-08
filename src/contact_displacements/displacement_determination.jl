using Colors

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

function stick_or_slip2(Nc, idxsΩ, w1_vec_full, w2_vec_full)
    m, n = Nc, length(idxsΩ)
    isstick = Matrix{Int}(undef, m, n)
    for i in 1:Nc
        for j in eachindex(idxsΩ)
            w1 = w1_vec_full[i,j]
            w2 = w2_vec_full[i,j]

            sub_w1 = w1_vec_full[i,j][end-2000:end]
            sub_w2 = w2_vec_full[i,j][end-2000:end]

            iszero1 = all(w1 .== 0.0)
            iszero2 = all(w2 .== 0.0)

            (min_w1, max_w1) = extrema(w1)
            (min_w2, max_w2) = extrema(w2)

            cond1 = !iszero1 && any(==(min_w1), sub_w1) && any(==(max_w1), sub_w1)
            cond2 = !iszero2 && any(==(min_w2), sub_w2) && any(==(max_w2), sub_w2)

            if iszero1 && iszero2
                isstick[i,j] = Int(0) # stick
            elseif cond1 || cond2 
                isstick[i,j] = Int(2) # totally slip (dissipates energy) 
            else 
                isstick[i,j] = Int(1) # partially slip (but not dissipates energy)
            end
        end
    end
    return isstick
end

