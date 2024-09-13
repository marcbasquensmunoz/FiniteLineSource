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
end
STSComputationContainers(n) = STSComputationContainers(zeros(n+1, n+1), zeros(n+1))

function initialize_containers(::SegmentToSegment, dps)
    N = map(dp -> dp.n, dps)
    unique_n = unique(N)
    map_n = [findall(x -> x==n, unique_n)[1] for n in N]
    (map(n -> STSComputationContainers(n), unique_n), map_n)
end

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp::DiscretizationParameters, containers::STSComputationContainers)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params
    @unpack P, M = containers

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

    for k in 0:n
        for s in 1:n+1
            P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    params = MeanSegToSegEvParams(setup)
    r_min, r_max = h_mean_lims(params) 
    guide(r, rb, params) = h_mean_sts(r*rb, params) * besselj(1/2, r) * imag(exp(im*r)) / r^(3/2)
    R̃, wz = adaptive_gk(r -> guide(r, rb, params), r_min/rb, r_max/rb)

    f(r̃, k, rb, m, c) = h_mean_sts(r̃*rb, params) * rb * besselj(k+1/2, m * r̃) * imag((im)^k * exp(im*c*r̃)) / r̃^(3/2)


    @inbounds for k in 0:n
        M[k+1] = dot(f.(R̃, k, rb, m, c), wz)
    end

    M .= P * M 
    return C .* M
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
