using LinearAlgebra, SparseArrays, Arpack

function craig_bampton_aplication(Mee, Mec, Mcc, Kee, Kec, Kcc, fe, fc; n_modes::Union{Int,Nothing}=nothing)

    #(Kee - ω² Mee) φ = 0
   # λs, Φ = eigs(Kee, Mee; nev=n_modes, which=:SR, ncv=3n_modes, maxiter=2000)  # Smallest Magnitude
    λs, Φ = eigs(Kee, Mee; nev=n_modes, which=:LR, sigma=1e-3, maxiter=3000)


    ω² = λs
    Φ = Φ

    #Ψ = -Kee⁻¹ Kec
    Ψ = -Kee \ Matrix(Kec)

    Max = Φ' * Mee * Ψ + Φ' * Mec
    Mxx = Ψ' * Mee * Ψ + Ψ' * Mec + Mec' * Ψ + Mcc
    Kxx = Kcc - Ψ' * Kee * Ψ

    fₐ = Φ' * fe
    fₓ = Ψ' * fe + fc

    #Ensamblado de matrices
    n_c = size(Kcc, 1)
    Id = I(n_modes)
    zerosmat = zeros(n_modes, n_c)

    M = [Id    Max;
         Max'  Mxx]

    K = [Diagonal(ω²)  zerosmat;
         zerosmat'         Kxx]

    #Matriz de transformacion T_cb
    
    T_CB = [Φ                  Ψ;
        zerosmat'  ones(n_c, n_c)]

    return M, Max, Mxx, K, Kxx, ω², fₐ, fₓ, Φ, Ψ, T_CB
end
