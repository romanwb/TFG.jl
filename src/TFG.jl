module TFG

#include("functions_4_test.jl")
#export my_f

using DifferentialEquations, GLMakie, Distributions, FFTW, LinearAlgebra, ContinuationSuite, Parameters

include("plots/themes.jl")
include("plots/plots.jl")

function __init__()
    apply_latex_style!()
end

include("duffing_params.jl")
include("fft.jl")
include("duffing_oscillator.jl")
# include("time_integration_comparison.jl") # (PROVISIONAL) Se eliminara tras separar en dos modulos
include("duffing_time_int.jl")
include("duffing_asymptotic.jl")
include("jacobian.jl")
include("friction_model.jl")
include("craig_bampton.jl")
include("ROM_data_processing.jl")
include("CB_cyclic.jl")


export DuffingParams, DuffingParamsContinuation, DuffingParamsTimeIntegration, DuffingParamsContinuationReal
export fft_matrices, system_matrix
export duffing, duffing_continuation, duffing_time_domain, duffing_time_domain_g, jacobian_matrix
# export DuffingTimeAsymp
export time_integration_values
export asymptotic_curve
export autodiff_jac, finite_diff_jac
export normal_force, tangencial_force, g_friction, g_friction_time
export craig_bampton_aplication
export read_rom_input, reduce_rom, save_rom_data, load_reduced_rom
export craig_bampton_aplication_sym

end