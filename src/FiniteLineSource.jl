module FiniteLineSource

export point_step_response
export precompute_parameters, compute_integral_throught_history!

using Parameters
using QuadGK
using LinearAlgebra
using FastGaussQuadrature
using LegendrePolynomials
using DSP
using Bessels
using StaticArrays
using SpecialFunctions
using Roots
using DataStructures
using FMM3D

include("nonhistory/interface.jl")
include("nonhistory/objects.jl")
include("nonhistory/bakhalov.jl")
include("nonhistory/point_to_point.jl")
include("nonhistory/segment_to_point.jl")
include("nonhistory/segment_to_segment.jl")
include("nonhistory/segment_to_segment_old.jl")
include("nonhistory/moving_point.jl")
include("nonhistory/moving_segment.jl")
include("nonhistory/moving_segment_to_segment.jl")
include("convolution.jl")
export PointToPoint, SegmentToPoint, SegmentToSegment, SegmentToSegmentOld, MovingPointToPoint, MovingSegmentToPoint, MovingSegmentToSegment
export convolve_step, step_response
export Constants
export precompute_parameters, compute_integral_throught_history!

include("approximations/mean_sts.jl")
include("approximations/point.jl")
include("approximations/mean_internal.jl")
export MeanSegToSegEvParams, PointEvalParams, InternalSegToSegEvParams

include("integration.jl")
export integrate

include("self_response.jl")
export compute_self_response

include("continuous/continuous.jl")
export compute_coefficients_through_history, precompute_matrices, legendre_coeffs

include("causal_non_history/asymptotic_integration.jl")
include("causal_non_history/block_method.jl")
include("causal_non_history/distance.jl")
include("causal_non_history/point_to_point.jl")
include("causal_non_history/line_to_point.jl")
include("causal_non_history/line_to_line.jl")
export AsymptoticContainers, evolve!, prepare_containers
end
