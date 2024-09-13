struct SegmentToSegment <: Setup
    D1
    H1
    D2
    H2
    σ
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

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp, P = nothing, M = nothing)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

    if isnothing(P)
        P = zeros(n+1, n+1)
    end
    if isnothing(M)
        M = zeros(n+1)
    end

    for k in 0:n
        for s in 1:n+1
            P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    params = MeanSegToSegEvParams(setup)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)
    guide(r) = h_mean_sts(r*rb) * besselj(1/2, r) * imag(exp(im*r)) / r^(3/2)
    R̃, wz = adaptive_gk(guide, r_min/rb, r_max/rb)

    f(r̃, k) = h_mean_sts(r̃*rb) * rb * besselj(k+1/2, m * r̃) * imag((im)^k * exp(im*c*r̃)) / r̃^(3/2)

    for k in 0:n
        M[k+1] = dot(f.(R̃, k), wz)
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
