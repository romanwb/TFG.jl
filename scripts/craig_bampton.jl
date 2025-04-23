using LinearAlgebra

function craig_bampton_aplication(Mee, Mec, Mce, Mcc, Kee, Kec, Kce, Kcc, fe, fc; n_modes::Int = 5)

    #(Kee - ω² Mee) φ = 0
    λs, Φ = eigen(Kee, Mee)
    idx = sortperm(λs)
    ω² = λs[idx][1:n_modes]
    Φ = Φ[:, idx[1:n_modes]]

    #Ψ = -Kee⁻¹ Kec
    Ψ = -Kee \ Kec

    #Definicion de terminos
    Max = Φ' * Mee * Ψ + Φ' * Mec
    Mxx = Ψ' * Mee * Ψ + Ψ' * Mec + Mce * Ψ + Mcc
    Kxx = Kcc - Ψ' * Kee * Ψ

    fa = Φ' * fe
    fx = Ψ' * fe + fc

    #Ensamblado de la matriz
    M = [Max'; Mxx]
    K = [Diagonal(ω²); Kxx]
    f = [fa; fx]

    #Matriz de transformacion T_cb
    n_c = size(Kcc, 1)
    T_CB = [
        Φ  Ψ;
        zeros(n_c, n_modes)  I(n_c)
    ]

    return M, K, f, Φ, Ψ, T_CB
end
