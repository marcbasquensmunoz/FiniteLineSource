import Bessels: besselj!
import .FiniteLineSource: DiscretizationParameters

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
    P::Matrix{T}
    M::Vector{T}
    aux::Vector{T}
    rule::Tuple{Vector{T}, Vector{T}, Vector{T}}
end
STSComputationContainers(n) = STSComputationContainers(zeros(n+1, n+1), zeros(n+1), zeros(n+1), QuadGK.kronrod(8))

function initialize_containers(::SegmentToSegment, dps)
    N = map(dp -> dp.n, dps)
    unique_n = unique(N)
    map_n = [findall(x -> x==n, unique_n)[1] for n in N]
    (map(n -> STSComputationContainers(n), unique_n), map_n)
end

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp::DiscretizationParameters, containers::STSComputationContainers)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params
    @unpack P, M, aux, rule = containers

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

    for k in 0:n
        for s in 1:n+1
            P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    sts_params = MeanSegToSegEvParams(setup)
    r_min, r_max = h_mean_lims(sts_params) 
    h_sts(r̃) = h_mean_sts(r̃*rb, sts_params)
    guide(r) = h_sts(r) * besselj(1/2, r) * imag(exp(im*r)) / r^(3/2)    
    R̃, wz = adaptive_gk(guide, r_min/rb, r_max/rb, rule=rule, rtol=1e-6)

    function f(r̃, N, rb, m, c, out)
        K = 0:N
        besselj!(out, K .+ 1/2, m * r̃)
        @. out *= h_sts(r̃) * rb * imag((im)^K * exp(im*c*r̃)) / r̃^(3/2)
    end

    J = zeros(n+1, length(R̃))

    for (i, r̃) in enumerate(R̃)
        f(r̃, n, rb, m, c, aux)
        @. J[:, i] = aux
    end

    for k in 0:n
        @views M[k+1] = dot(J[k+1, :], wz)
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

