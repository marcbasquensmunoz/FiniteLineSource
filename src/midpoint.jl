
@with_kw struct MidPointParams
    D
    H
    z
    σ
    rb

    rmin = σ/rb
    rB = sqrt(σ^2 + (z - D - H)^2 )/rb 
    rT = sqrt(σ^2 + (z - D)^2     )/rb 

    r1 = (D < z && z < D+H) ? rmin : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)
end

function h(r, σ, rb, r1, r2, r3)
    if r < r1
        0
    elseif r < r2
        2r*rb^2/sqrt((r*rb)^2-σ^2)
    elseif r < r3
        r*rb^2/sqrt((r*rb)^2-σ^2)
    else
        0
    end
end

function midpoint_evaluation(params::MidPointParams)
    @unpack D, H, z, σ, rb, r1, r2, r3 = params

    h_midpoint(r) = h(r, σ, rb, r1, r2, r3)
    r = r1:0.01:r3

    return h_midpoint, r
end
