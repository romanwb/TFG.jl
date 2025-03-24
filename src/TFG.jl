module TFG

#include("functions_4_test.jl")
#export my_f

using DifferentialEquations, GLMakie, Distributions, FFTW, LinearAlgebra, ContinuationSuite, Parameters
include("duffing_params.jl")
include("fft.jl")
include("duffing_oscillator.jl")
# include("time_integration_comparison.jl") # (PROVISIONAL) Se eliminara tras separar en dos modulos
include("duffing_time_int.jl")
include("duffing_asymptotic.jl")



export DuffingParams, DuffingParamsContinuation, DuffingParamsTimeIntegration
export fft_matrices, system_matrix
export duffing, duffing_continuation, duffing_time_domain, jacobian_matrix
# export DuffingTimeAsymp
export time_integration_values
export asymptotic_curve


end