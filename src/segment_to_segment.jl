@with_kw struct SegmentToSegment{T <: Number} <: Setup @deftype T
    D1
    H1
    D2
    H2
    r
end

function Preallocation(::SegmentToSegment, params::Constants) 
    @unpack segment_points, line_points = params
    P = [zeros(i+1, i+1) for i in segment_points]
    M = [zeros(i+1) for i in segment_points]
    Preallocation(P=P, R=[zeros(0, 0) for i in segment_points], M=M)
end

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp, P, R, M)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

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
    @unpack kg, α = params
    @unpack D1, H1, D2, H2, r = setup
    I(d) = sqrt(r^2+d^2) + d * log(sqrt(r^2+d^2) - d)
    1/(4π*kg*H2) * (I(D1+H1-D2-H2) + I(D1-D2) - I(D1-H2-D2) - I(D1+H1-D2))
end

function analytical_test(setup::SegmentToSegment; params::Constants, t) 
    @unpack D1, H1, D2, H2, r = setup
    @unpack kg, α = params
    quadgk(z -> quadgk(ζ -> 1/(4π*kg*H2*sqrt(r^2+(z-ζ)^2)) * erfc(sqrt(r^2+(z-ζ)^2)/sqrt(4α*t)), D1, H1+D1, atol=10^-16, order=20)[1], D2, D2+H2, atol=10^-16, order=20)
end

function has_heatwave_arrived(setup::SegmentToSegment; params::Constants, t)
    @unpack D1, H1, D2, H2, r = setup
    @unpack α = params
    threshold = 8

    if D2 >= D1 && D2 <= D1+H1 || D1 >= D2 && D1 <= D2+H2
        d2 = r^2
    else 
        if D2 > D1+H1
            d2 = r^2 + (D2-D1-H1)^2 
        else 
            d2 = r^2 + (D1-D2-H2)^2 
        end
    end

    d2 / (2α*t) < threshold^2
end

