using LinearAlgebra

function craig_bampton_projection(
    Mee, Mec, Mce, Mcc,
    Kee, Kec, Kce, Kcc,
    fe, fc;
    n_modes::Int = 5
)

    # 1. Autovalores del sistema full-stuck: (Kee - ω² Mee) φ = 0
    λs, Φ = eigen(Kee, Mee)
    idx = sortperm(λs)
    ω² = λs[idx][1:n_modes]
    Φ = Φ[:, idx[1:n_modes]]

    # 2. Modos estáticos Ψ = -Kee⁻¹ Kec
    Ψ = -Kee \ Kec

    # 3. Proyecciones de masa, rigidez y fuerzas
    Max = Φ' * Mee * Ψ + Φ' * Mec
    Mxx = Ψ' * Mee * Ψ + Ψ' * Mec + Mce * Ψ + Mcc
    Kxx = Kcc - Ψ' * Kee * Ψ

    fa = Φ' * fe
    fx = Ψ' * fe + fc

    # 4. Ensamblado de matrices proyectadas
    M = [Max'; Mxx]
    K = [Diagonal(ω²); Kxx]
    f = [fa; fx]

    # 5. Matriz de transformación T_CB
    n_c = size(Kcc, 1)
    T_CB = [
        Φ  Ψ;
        zeros(n_c, n_modes)  I(n_c)
    ]

    return M, K, f, Φ, Ψ, T_CB
end
