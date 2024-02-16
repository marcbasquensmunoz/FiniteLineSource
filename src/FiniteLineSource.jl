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

include("BoreholeSegment.jl")
include("segment_response.jl")
include("segment_T_field.jl")

include("iterative.jl")
include("point_to_point.jl")
include("segment_to_point.jl")
include("segment_to_segment.jl")
export PointToPoint, SegmentToPoint, SegmentToSegment
export Constants, Preallocation
export precompute_parameters, compute_integral_throught_history!

end
