
@with_kw struct MeanSegToSegEvParams{T <: Number} @deftype T
    D1
    H1
    D2
    H2
    σ

    rmin = σ

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin
end
MeanSegToSegEvParams(p::SegmentToSegment) = MeanSegToSegEvParams(D1=p.D1, H1=p.H1, D2=p.D2, H2=p.H2, σ=p.σ)
MeanSegToSegEvParams(p::MovingSegmentToSegment) = MeanSegToSegEvParams(D1=p.D1, H1=p.H1, D2=p.D2, H2=p.H2, σ=p.σ)
transpose(p::MeanSegToSegEvParams) = MeanSegToSegEvParams(D1=p.D2, H1=p.H2, D2=p.D1, H2=p.H1, σ=p.σ)

function L(r, params::MeanSegToSegEvParams)
    @unpack D1, H1, D2, H2, σ, r1, r2, r3, r4 = params

    if r <= r1
        return 0.
    elseif r < r2
        return r * (D1 - D2 + H1) / sqrt(r^2-σ^2) + r
    elseif r < r3
        return  r * min(H1, H2) / sqrt(r^2-σ^2)
    elseif r < r4
        return r * (D2 - D1 + H2) / sqrt(r^2-σ^2) - r
    else
        return 0.
    end
end

h_mean_sts(r, p::MeanSegToSegEvParams) = (L(r, p) + L(r, transpose(p))) / p.H2
function h_mean_lims(params)
    paramsT = transpose(params)
    r_min = Base.max(params.r1, paramsT.r1)
    r_max = Base.max(params.r4, paramsT.r4)
    return r_min, r_max
end
