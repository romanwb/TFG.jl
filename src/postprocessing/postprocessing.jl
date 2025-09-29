struct ProjectorMaker{T}
    E::Matrix{T}
    m::Int
end

function ProjectorMaker(N::Integer, H::Integer)
    E, Eᴴ = HBMContinuation.dft_matrices(N, H)
    m = 2H+1
    return ProjectorMaker(E, m)
end

function make_FRF(op::ProjectorMaker, Xⱼ, Aⱼ; dim::Int=3)
    m = op.m
    Nc = size(Xⱼ, 1) ÷ (dim * m)
    Na = size(Aⱼ, 1) ÷ m

    x = Vector{Matrix{Float64}}(undef, 3Nc)
    x_max = Vector{Vector{Float64}}(undef, dim * Nc)

    a = Vector{Matrix{Float64}}(undef, Na)
    a_max = Vector{Vector{Float64}}(undef, Na)

    for i in 1:3Nc
        r1 = (m-1)*(i-1) + i
        r2 = r1 + m - 1

        x[i] = op.E * Xⱼ[r1:r2, :]
        x_max[i] = vec(maximum(abs.(x[i]), dims = 1))
    end

    for i in 1:Na
        r1 = (m-1)*(i-1) + i
        r2 = r1 + m - 1

        a[i] = op.E * Aⱼ[r1:r2, :]
        a_max[i] = vec(maximum(abs.(a[i]), dims = 1))
    end

    return x_max, a_max, x, a
end

