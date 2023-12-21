module FiniteLineSource

export BoreholeSegment, segment_to_segment_step_response, T_ls
export point_step_response
export precompute_parameters, fevolve_through_history!, compute_integral_and_advance_one_step!, compute_integral_throught_history!

using SpecialFunctions
using QuadGK

include("BoreholeSegment.jl")
include("segment_response.jl")
include("segment_T_field.jl")
include("bakhvalov.jl")

end
