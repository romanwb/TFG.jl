using LinearAlgebra, SparseArrays, Arpack

function craig_bampton_aplication(Mee, Mec, Mcc, Kee, Kec, Kcc, fe, fc; n_modes::Union{Int,Nothing}=nothing)

    #(Kee - ω² Mee) φ = 0
    λs, Φ = eigs(Kee, Mee; nev=n_modes, which=:SM)  # Smallest Magnitude
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
    Id = ones(n_modes, n_modes)
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

function craig_bampton(Mee, Mec, Mcc, Kee, Kec, Kcc, fe, fc; n_modes::Union{Int,Nothing}=nothing)

       #(Kee - ω² Mee) φ = 0
       λs, Φ_full = eigen(Kee, Mee)
       idx = sortperm(λs) #se obtienen posiciones de los autovalores ordenados
   
       n_total = length(λs)
       n_modes = isnothing(n_modes) ? n_total : min(n_modes, n_total)
   
       ω² = λs[idx][1:n_modes]
       Φ = Φ_full[:, idx[1:n_modes]]
   
       #Ψ = -Kee⁻¹ Kec
       Ψ = -Kee \ Kec
   
       Max = Φ' * Mee * Ψ + Φ' * Mec
       Mxx = Ψ' * Mee * Ψ + Ψ' * Mec + Mec' * Ψ + Mcc
       Kxx = Kcc - Ψ' * Kee * Ψ
   
       fₐ = Φ' * fe
       fₓ = Ψ' * fe + fc
   
       #Ensamblado de matrices
       n_c = size(Kcc, 1)
       Id = spdiagm(0 => ones(n_modes))
       zerosmat = spzeros(n_modes, n_c)
   
       M = [sparse(Id)     sparse(Max);
            sparse(Max')   sparse(Mxx)]
   
       K = [spdiagm(0 => ω²)  zerosmat;
            zerosmat'         sparse(Kxx)]
   
       #Matriz de transformacion T_cb
       
       T_CB = [Φ                  Ψ;
           spzeros(n_c, n_modes)  spdiagm(0 => ones(n_c))]
   
       return M, Max, Mxx, K, Kxx, ω², fₐ, fₓ, Φ, Ψ, T_CB
   end 


# Test 
m=1.0
k=10.0

Mee = [m 0.0;
       0.0 m]

Mec = [0.0 0.0 0.0;
       0.0 0.0 0.0]

Mcc = [m 0.0 0.0;
       0.0 m 0.0;
       0.0 0.0 m]

Kee = [k -k;
       -k 2k]

Kec = [0.0 0.0 0.0;
       -k 0.0 0.0]


Kcc = [2k -k 0.0;
       -k 2k -k;
       0  -k  2k]

fe = [1.0; 1.0]

fc = [0.0; 0.0; 0.0]

n_modes = 2

M, Max, Mxx, K, Kxx, ω², fₐ, fₓ, Φ, Ψ, T_CB = craig_bampton(Mee, Mec, Mcc, Kee, Kec, Kcc, fe, fc; n_modes)
Ψ


fₓ
# Test SparseArrays

Mee = sparse([m 0.0;
              0.0 m])

Mec = spzeros(2, 3)

Mcc = sparse([m 0.0 0.0;
              0.0 m 0.0;
              0.0 0.0 m])

Kee = sparse([k -k;
              -k 2k])

Kec = sparse([0.0 0.0 0.0;
              -k 0.0 0.0])

Kcc = sparse([2k -k 0.0;
              -k 2k -k;
              0.0 -k 2k])

fe = [1.0; 1.0]

fc = [0.0; 0.0; 0.0]

n_modes = 1

M, Max, Mxx, K, Kxx, ω², fₐ, fₓ, Φ, Ψ, T_CB = craig_bampton_aplication(Mee, Mec, Mcc, Kee, Kec, Kcc, fe, fc; n_modes)

Kee_inv= inv(Kee)
Ψ