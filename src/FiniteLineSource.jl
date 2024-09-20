module FiniteLineSource

export point_step_response
export precompute_parameters, fevolve_through_history!, compute_integral_and_advance_one_step!, compute_integral_throught_history!

using Parameters
using QuadGK
using LinearAlgebra
using Cubature
using FastGaussQuadrature
using LegendrePolynomials
using DSP
using Bessels
using StaticArrays

include("nonhistory/interface.jl")
include("nonhistory/objects.jl")
include("nonhistory/bakhalov.jl")
include("nonhistory/point_to_point.jl")
include("nonhistory/segment_to_point.jl")
include("nonhistory/segment_to_segment.jl")
include("nonhistory/segment_to_segment_old.jl")
include("nonhistory/moving_point.jl")
include("nonhistory/moving_segment.jl")
include("convolution.jl")
export PointToPoint, SegmentToPoint, SegmentToSegment, SegmentToSegmentOld, MovingPointToPoint, MovingSegmentToPoint
export convolve_step, step_response
export Constants, Preallocation
export precompute_parameters, compute_integral_throught_history!

include("approximations/mean_sts.jl")
include("approximations/midpoint.jl")
include("approximations/mean_internal.jl")
export mean_sts_evaluation, midpoint_evaluation, mean_internal_evaluation
export MeanSegToSegEvParams, MidPointParams, InternalSegToSegEvParams

include("integration.jl")
export integrate
export IntegrationSegment

include("self_response.jl")
export compute_self_response

include("continuous/continuous.jl")
export compute_coefficients_through_history, precompute_matrices, legendre_coeffs

end
