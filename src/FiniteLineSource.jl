module FiniteLineSource

export BoreholeSegment, segment_to_segment_step_response, T_ls
export compute_integral, point_step_response
export discretization_parameters, frequency_parameters, frequency_parameters_2, evolve!

using SpecialFunctions
using QuadGK

include("BoreholeSegment.jl")
include("segment_response.jl")
include("segment_T_field.jl")
include("bakhvalov.jl")

end
