using HBMContinuation, GLMakie, LinearAlgebra, Julianim, TFG, JLD2, SparseArrays
use_formal_theme!()
@load "scripts/realistic_models/turbine_blade_model/plot_data/EO1/7x7_ff1.jld2" Ω Xⱼ Aⱼ p
    H = p.H
    Nx = Int(length(Xⱼ[:,1])/(11))
    Nc = Int(Nx/3)
    n = Int(sqrt(Nc))
    N = length(Xⱼ[1,:])
    X₁ = Xⱼ[1:(2H + 1), :]
    X=zeros(Nx,2H+1,N)
    for i in 1:Nx
        X[i,:,:] = Xⱼ[(i-1)*(2H + 1)+1:i*(2H + 1), :]
    end

    E, Eᴴ = HBMContinuation.dft_matrices(2^10, H)
    physic_x⃗ = zeros(Nx,N)
    for i in 1:Nx
    physic_x⃗[i,:] = vec(maximum(abs.(E * X[i,:,:]), dims = 1))
    end
    lines(Ω, physic_x⃗[1,:])
    
    physic_x⃗1 = physic_x⃗[:,1]
    physic_x⃗1x = zeros(Nc)
    physic_x⃗1y = zeros(Nc)
    physic_x⃗1z = zeros(Nc)
    for i in 1:Nc
    physic_x⃗1x[i] = physic_x⃗1[3(i-1)+1]
    physic_x⃗1y[i] = physic_x⃗1[3(i-1)+2]
    physic_x⃗1z[i] = physic_x⃗1[3(i-1)+3]
    end

    Xcoord = reshape(physic_x⃗1x, n, n)
    Ycoord = reshape(physic_x⃗1y, n, n)
    Zcoord = reshape(physic_x⃗1z, n, n)
    


function plot_deformed_grid3d(Ux, Uy, Uz; ny::Int=5, nx::Int=5,
                              X0=nothing, Y0=nothing, Z0=nothing, scale::Real=1.0)

    UxM = Ux isa AbstractMatrix ? Ux : reshape(Ux, ny, nx)
    UyM = Uy isa AbstractMatrix ? Uy : reshape(Uy, ny, nx)
    UzM = Uz isa AbstractMatrix ? Uz : reshape(Uz, ny, nx)

    if X0 === nothing || Y0 === nothing || Z0 === nothing
        xv = range(0, 1, length=nx)
        yv = range(0, 1, length=ny)
        X0M = [x for y in yv, x in xv]     # ny×nx
        Y0M = [y for y in yv, x in xv]
        Z0M = zeros(ny, nx)
    else
        X0M = X0 isa AbstractMatrix ? X0 : reshape(X0, ny, nx)
        Y0M = Y0 isa AbstractMatrix ? Y0 : reshape(Y0, ny, nx)
        Z0M = Z0 isa AbstractMatrix ? Z0 : reshape(Z0, ny, nx)
    end

    X = X0M .+ scale .* UxM
    Y = Y0M .+ scale .* UyM
    Z = Z0M .+ scale .* UzM

    fig = Figure()
    ax  = Axis3(fig[1, 1], title = "Malla 3D: base vs deformada")

    for j in 1:nx
        lines!(ax, X0M[:, j], Y0M[:, j], Z0M[:, j], color=:gray, linewidth=1)
    end
    for i in 1:ny
        lines!(ax, X0M[i, :], Y0M[i, :], Z0M[i, :], color=:gray, linewidth=1)
    end

    for j in 1:nx
        lines!(ax, X[:, j], Y[:, j], Z[:, j], color=:blue, linewidth=2)
    end
    for i in 1:ny
        lines!(ax, X[i, :], Y[i, :], Z[i, :], color=:blue, linewidth=2)
    end
    scatter!(ax, vec(X), vec(Y), vec(Z), color=:blue, markersize=8)

    fig
end


ny, nx = n, n


fig = plot_deformed_grid3d(Xcoord, Ycoord, Zcoord; ny=ny, nx=nx, scale=1.0)
display(fig)
