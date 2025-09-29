module TFG

#include("functions_4_test.jl")
#export my_f

using DifferentialEquations, GLMakie, Distributions, FFTW, LinearAlgebra, ContinuationSuite, Parameters

include("plots/constants_theme.jl")
include("plots/plot_theme.jl")

include("duffing_params.jl")
include("fft.jl")
include("duffing_oscillator.jl")

include("duffing_time_int.jl")

include("jacobian.jl")
include("friction_model.jl")
include("craig_bampton.jl")
include("ROM_data_processing.jl")

include("mat_jld2.jl")

include("postprocessing/postprocessing.jl")

include("contact_displacements/contact_displacements.jl")


export CircularContainer, COLORS, GREYS, MARKERS
export use_formal_theme!

export DuffingParams, DuffingParamsContinuation, DuffingParamsTimeIntegration, DuffingParamsContinuationReal
export fft_matrices, system_matrix
export duffing, duffing_continuation, duffing_time_domain, duffing_time_domain_g, jacobian_matrix

export time_integration_values

export autodiff_jac, finite_diff_jac
export normal_force, tangencial_force, g_friction, g_friction_time
export craig_bampton_aplication, craig_bampton_aplication_sym
export read_rom_input, reduce_rom, save_rom_data, load_reduced_rom

export convert_mat_to_jld2

export ProjectorMaker, make_FRF

export study_1contact, study_fullplane, make_indexs, full_mesh, stick_or_slip, plot_slip_stick

end