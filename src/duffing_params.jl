"""
    duffing_params

Estructura que encapsula los parámetros necesarios para resolver el problema del oscilador
de Duffing usando la transformada de Fourier.

- 'N': puntos en los que se evalua la funcion respectivamente.
- 'H': numero de armónicos retenidos en la serie de Fourier.
- 'E', 'Eᴴ': matrices para la transformada directa e inversa (creadas con 'fft_matrices').
- 'x₀': vector de condiciones iniciales de dimensión 2H + 1.
- 'f̂': el vector de Fourier para la fuerza (también de dimensión '2H + 1').
- 'g', 'dg': funciones no lineales (por ejemplo, 'g(x)=x^3', 'dg(x)=3x^2').
- 'Ω', 'ξ', 'ϵ': parámetros físicos del problema.

Se valida que los vectores 'x₀' y 'f̂' tengan siempre longitud '2H+1'.
"""

using Parameters
@with_kw struct DuffingParams

    N::Int
    H::Int

    Ω::Float64
    ξ::Float64
    ϵ::Float64

    E::Matrix{Float64}
    Eᴴ::Matrix{Float64}

    x₀::Vector{Float64}

    @assert length(x₀) == 2H + 1 "x₀ debe tener dimensión 2H+1, pero tiene $(length(x₀))"

    f̂::Vector{Float64}

    @assert length(f̂) == 2H + 1 "f̂ debe tener dimensión 2H+1, pero tiene $(length(f̂))"

    g::Function
    dg::Function

end


using Parameters
@with_kw struct DuffingParamsContinuation

    N::Int
    H::Int
    
    ξ::Float64
    ϵ::Float64

    E::Matrix{Float64}
    Eᴴ::Matrix{Float64}

    x₀::Vector{Float64}

    @assert length(x₀) == 2H + 1 "x₀ debe tener dimensión 2H+1, pero tiene $(length(x₀))"

    f̂::Vector{Float64}

    @assert length(f̂) == 2H + 1 "f̂ debe tener dimensión 2H+1, pero tiene $(length(f̂))"

    g::Function
    dg::Function

end