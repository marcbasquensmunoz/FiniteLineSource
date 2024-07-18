
@with_kw struct MeanSegToSegEvParams{T <: Number} @deftype T
    D1
    H1
    D2
    H2
    σ

    rmin = σ

    rLRp = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLLp = sqrt(σ^2 + (D2 - D1)^2          )
    rULp = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rURp = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    rLR = H1 > 0 ? (H2 > 0 ? rLRp : rURp) : (H2 > 0 ? rLLp : rULp)
    rLL = H1 > 0 ? (H2 > 0 ? rLLp : rULp) : (H2 > 0 ? rLRp : rURp)
    rUL = H1 > 0 ? (H2 > 0 ? rULp : rLLp) : (H2 > 0 ? rURp : rLRp)
    rUR = H1 > 0 ? (H2 > 0 ? rURp : rLRp) : (H2 > 0 ? rULp : rLLp)

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin
end
MeanSegToSegEvParams(p::SegmentToSegment) = MeanSegToSegEvParams(D1=p.D1, H1=p.H1, D2=p.D2, H2=p.H2, σ=p.σ)

function swap!(a, b)
    aux = b
    b = a
    a = aux
end

transpose(p::MeanSegToSegEvParams) = MeanSegToSegEvParams(D1=p.D2, H1=p.H2, D2=p.D1, H2=p.H1, σ=p.σ)

function L(r, params::MeanSegToSegEvParams)
    @unpack D1, H1, D2, H2, σ, r1, r2, r3, r4 = params

    if r <= r1
        0.
    elseif r < r2
        r * (D1 - D2 + H1) / sqrt(r^2-σ^2) + r
    elseif r < r3
        r * min(H1, H2) / sqrt(r^2-σ^2)
    elseif r < r4
        r * (D2 - D1 + H2) / sqrt(r^2-σ^2) - r
    else
        0.
    end
end

function mean_sts_evaluation(params::MeanSegToSegEvParams)
    paramsT = transpose(params)
    f(r) = L(r, params) + L(r, paramsT)
    h_mean_sts(r) = f(r) / params.H2
    r_min = Base.max(params.r1, paramsT.r1)
    r_max = Base.max(params.r4, paramsT.r4)
    return h_mean_sts, r_min, r_max
end

function compute_double_integral_sts(s::SegmentToSegment, f::Function)
    params = MeanSegToSegEvParams(s)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)
    x, w = adaptive_gk(f, r_min, r_max)
    dot(f.(x), w)
end