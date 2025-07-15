using LinearAlgebra, SparseArrays, Arpack, HBMContinuation, JLD2

     function f_k_from_wave(fe::AbstractVector, k::Int, N::Int, EO::Int)
          if k % N == EO % N
          return ComplexF64.(fe)
          else
          return zeros(ComplexF64, size(fe))
          end
     end

     # Gram-Schmidt
     function orthonormalize_M(Φ::Matrix{ComplexF64}, M::AbstractMatrix{<:Complex})
          n = size(Φ, 2)
          Φ̃ = similar(Φ)
      
          for j in 1:n
              v = Φ[:, j]
      
              for i in 1:j-1
                  proj = (Φ̃[:, i]' * M * v)
                  v -= Φ̃[:, i] * proj
              end
      
              normM = sqrt(real(v' * M * v))
              Φ̃[:, j] = v / normM
          end
      
          return Φ̃
      end      
      

   function craig_bampton_aplication_sym(Mee, Mec, Mcc, MCee, Kee, Kec, Kcc, KCee, fe, fc, Rx; k::Int, N::Int, EO::Int, n_modes::Union{Int,Nothing}=nothing)
 
     Ne = size(Mee, 1)
     Nc = size(Mcc, 1)
 
     Zec = zeros(Ne, Nc)
     Zcc = zeros(Nc, Nc)
     Mc = zeros(Ne+Nc, Ne+Nc)
     Mc[1:size(MCee,1),1:size(MCee,2)] = MCee

     Kc = zeros(Ne+Nc, Ne+Nc)
     Kc[1:size(KCee,1),1:size(KCee,2)] = KCee
 
     M = [Mee  Mec;
          Mec' Mcc]
 
     K = [Kee  Kec;
          Kec' Kcc]
 
     θ = 2π * k / N
     phase = exp(im * θ)
     Mᵏ = transpose(Mc) * conj(phase) + M + Mc * phase
     Kᵏ = transpose(Kc) * conj(phase) + K + Kc * phase
 
     Keeᵏ = Kᵏ[1:Ne, 1:Ne]
     Kecᵏ = Kᵏ[1:Ne, Ne+1:end]
     Kceᵏ = Kecᵏ'
     Kccᵏ = Kᵏ[Ne+1:end, Ne+1:end]
 
     Meeᵏ = Mᵏ[1:Ne, 1:Ne]
     Mecᵏ = Mᵏ[1:Ne, Ne+1:end]
     Mceᵏ = Mecᵏ'
     Mccᵏ = Mᵏ[Ne+1:end, Ne+1:end]
 
     λs, Φ = eigs(Keeᵏ, Meeᵏ; nev=n_modes, which=:LR, sigma=0.0, maxiter=200000)
 
     ωᵏ² = λs
     Φᵏ = Φ
     Φᵏ = orthonormalize_M(Φᵏ, Meeᵏ)


     Ψᵏ = -Keeᵏ \ Matrix(Kecᵏ)
 
     T_CB = [Φᵏ                      Ψᵏ;
             zeros(Nc, size(Φᵏ,2))   I(Nc)]
 
     Φᵏᴴ = Φᵏ'
     Ψᵏᴴ = Ψᵏ'
 
     Maxᵏ = Φᵏᴴ * Meeᵏ * Ψᵏ + Φᵏᴴ * Mecᵏ
     Mxxᵏ = Ψᵏᴴ * Meeᵏ * Ψᵏ + Ψᵏᴴ * Mecᵏ + Mceᵏ * Ψᵏ + Mccᵏ
     Kxxᵏ = Kccᵏ - Ψᵏᴴ * Keeᵏ * Ψᵏ

     fe_full = zeros(Ne)
     # .rowval -> vector con indices de las filas de elementos no nulos
     # .nzval -> vector con valores no nulos
     fe_full[fe.rowval[1]] = fe.nzval[1]

     fᵏ = f_k_from_wave(fe_full, k, N, EO)
     faᵏ = Φᵏᴴ * fᵏ
     fxᵏ = Ψᵏᴴ * fᵏ

     return M, Maxᵏ, Mxxᵏ, K, Kxxᵏ, real(ωᵏ²), faᵏ, fxᵏ, Φᵏ, Ψᵏ, T_CB
 end
