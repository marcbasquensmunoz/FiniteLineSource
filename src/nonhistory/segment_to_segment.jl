import Bessels: besselj!
import .FiniteLineSource: DiscretizationParameters
import .Bessels: besselj

struct SegmentToSegment{T <: Number} <: Setup
    D1::T
    H1::T
    D2::T
    H2::T
    σ::T
end
function SegmentToSegment(;D1, H1, D2, H2, σ)
    lowest = min(min(D1, D1+H1), min(D2, D2+H2))
    if lowest < 0
        lowest = abs(lowest)
        nD1, nH1 = min(D1+lowest, D1+H1+lowest), max(D1+lowest, D1+H1+lowest) - min(D1+lowest, D1+H1+lowest)
        nD2, nH2 = min(D2+lowest, D2+H2+lowest), max(D2+lowest, D2+H2+lowest) - min(D2+lowest, D2+H2+lowest)
        SegmentToSegment(nD1, nH1, nD2, nH2, σ)
    else 
        return SegmentToSegment(D1, H1, D2, H2, σ)
    end
end

struct STSComputationContainers{T <: Number} <: ComputationContainers
    J::Matrix{T}
    P::Matrix{T}
    M::Vector{T}
    aux::Vector{T}
    segment_buffer::Union{Vector{QuadGK.Segment{T, T, T}}, Nothing}
end
STSComputationContainers(n) = STSComputationContainers(zeros(n+1, 2*n^3), zeros(n+1, n+1), zeros(n+1), zeros(n+1), nothing)

function initialize_containers(::SegmentToSegment, dps)
    N = map(dp -> dp.n, dps)
    unique_n = unique(N)
    map_n = [findfirst(x -> x==n, unique_n) for n in N]
    (map_n, map(n -> STSComputationContainers(n), unique_n))
end

function initialize_buffer(setup::SegmentToSegment, rb)
    sts_params = MeanSegToSegEvParams(setup)
    r_min, r_max = h_mean_lims(sts_params) 
    h_sts(r̃) = h_mean_sts(r̃*rb, sts_params)
    guide(r) = h_sts(r) * besselj(1/2, r) * imag(exp(im*r)) / r^(3/2)  
    _,_, buffer = quadgk_segbuf(guide, r_min/rb, r_max/rb, order = 20)
    return buffer
end

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp::DiscretizationParameters, containers::STSComputationContainers, buffer=nothing, rtol=sqrt(eps()), atol=0.)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params
    @unpack J, P, M, aux = containers

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

    for k in 0:n
        for s in 1:n+1
            @inbounds P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    sts_params = MeanSegToSegEvParams(setup)
    #sts_params_T = transpose(sts_params)

    r_min, r_max = h_mean_lims(sts_params) 
    #h_sts(r̃) = (L(r̃*rb, sts_params) + L(r̃*rb, sts_params_T)) / sts_params.H2 
    h_sts(r̃) = h_mean_sts(r̃*rb, sts_params)
    guide(r̃) = h_sts(r̃) * besselj(1/2, r̃) * imag(exp(im*r̃)) / r̃^(3/2)  
    R̃, wz = adaptive_nodes_and_weights(guide, r_min/rb, r_max/rb, n = 20, buffer = buffer, rtol=rtol, atol=atol)
    
    function f(r̃, N, rb, m, c, out)
        besselj!(out, 1/2:(N+1/2), m * r̃)
        eval = h_sts(r̃) * rb * exp(im*c*r̃) / r̃^(3/2)
        value = (imag(eval), imag(im*eval), -imag(eval), -imag(im*eval))
        for k in 0:N
            @inbounds out[k+1] *= value[(k % 4) + 1]
        end
    end
    
    for (i, r̃) in enumerate(R̃)
        f(r̃, n, rb, m, c, aux)
        for k in 1:n+1
            @inbounds J[k, i] = aux[k]
        end
    end

    for k in 0:n
        @inbounds @views M[k+1] = dot(J[k+1, 1:length(wz)], wz)
    end

    mul!(aux, P, M)
    return C .* aux
end

function constant_integral(setup::SegmentToSegment; params::Constants)
    @unpack kg = params
    @unpack D1, H1, D2, H2, σ = setup
    β(d) = sqrt(σ^2+d^2) + d * log(sqrt(σ^2+d^2) - d)
    1/(4π*kg*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1-H2-D2) - β(D1+H1-D2))
end

function has_heatwave_arrived(setup::SegmentToSegment; params::Constants, t)
    @unpack D1, H1, D2, H2, σ = setup
    @unpack α = params
    threshold = 8

    if D2 >= D1 && D2 <= D1+H1 || D1 >= D2 && D1 <= D2+H2
        d2 = σ^2
    else 
        if D2 > D1+H1
            d2 = σ^2 + (D2-D1-H1)^2 
        else 
            d2 = σ^2 + (D1-D2-H2)^2 
        end
    end

    d2 / (2α*t) < threshold^2
end

