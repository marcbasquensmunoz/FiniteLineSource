
@with_kw struct PointEvalParams
    D
    H
    z
    σ

    rmin = σ
    rB = sqrt(σ^2 + (z - D - H)^2 )
    rT = sqrt(σ^2 + (z - D)^2     )

    r1 = (D < z && z < D+H) ? rmin : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)
end

PointEvalParams(s::SegmentToPoint) = PointEvalParams(D=s.D, H=s.H, z=s.z, σ=s.σ)

function h_point(r, params::PointEvalParams)
    @unpack σ, r1, r2, r3 = params
    if r <= r1
        return 0.
    elseif r < r2
        return 2r/sqrt(r^2-σ^2)
    elseif r < r3
        return r/sqrt(r^2-σ^2)
    else
        return 0.
    end
end

function h_point_lims(params::PointEvalParams)
    @unpack r1, r3 = params
    return r1, r3
end
