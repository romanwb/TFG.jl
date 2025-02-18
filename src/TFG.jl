module TFG

include("functions_4_test.jl")
export my_f

using DifferentialEquations, CairoMakie, Distributions
include("time_integration_w_cte.jl")
include("time_int_comparison.jl")

end
