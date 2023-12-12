module FiniteLineSource

export BoreholeSegment, segment_to_segment_step_response, T_ls, compute_integral

include("BoreholeSegment.jl")
include("segment_response.jl")
include("segment_T_field.jl")
include("bakhvalov.jl")

end
