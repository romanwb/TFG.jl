module TFG

#include("functions_4_test.jl")
#export my_f

using DifferentialEquations, GLMakie, Distributions, FFTW, LinearAlgebra, ContinuationSuite, Parameters
include("duffing_params.jl")
include("fft.jl")
include("duffing_oscillator.jl")
include("time_integration_comparison.jl")



export DuffingParams, DuffingParamsContinuation
export fft_matrices, system_matrix
export duffing, duffing_continuation, jacobian_matrix
export DuffingTimeAsymp



end