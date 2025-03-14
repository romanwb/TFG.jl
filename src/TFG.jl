module TFG

#include("functions_4_test.jl")
#export my_f

using DifferentialEquations, GLMakie, Distributions, FFTW, LinearAlgebra, ContinuationSuite, Parameters
include("duffing_params.jl")
include("fft.jl")
include("duffing_oscillator.jl")



export DuffingParams
export fft_matrices, system_matrix
export duffing, jacobian_matrix




end