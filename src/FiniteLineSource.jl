module FiniteLineSource

export BoreholeSegment, segment_to_segment_step_response, T_ls
export point_step_response
export precompute_parameters, fevolve_through_history!, compute_integral_and_advance_one_step!, compute_integral_throught_history!

using Parameters
using SpecialFunctions
using QuadGK
using LinearAlgebra
using Cubature
using FastGaussQuadrature
using LegendrePolynomials
using DSP

include("BoreholeSegment.jl")
include("segment_response.jl")
include("segment_T_field.jl")


include("bakhalov.jl")
include("point_to_point.jl")
include("segment_to_point.jl")
include("segment_to_segment.jl")
include("moving_point.jl")
include("convolution.jl")
export PointToPoint, SegmentToPoint, SegmentToSegment, MovingPointToPoint
export convolve_step, step_response
export Constants, Preallocation
export precompute_parameters, compute_integral_throught_history!

include("mean_sts.jl")
include("midpoint.jl")
export mean_sts_evaluation, midpoint_evaluation
export MeanSegToSegEvParams, MidPointParams

include("integration.jl")
export integrate

include("self_response.jl")
export compute_self_response

end
