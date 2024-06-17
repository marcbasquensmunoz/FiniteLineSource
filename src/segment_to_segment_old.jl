@with_kw struct SegmentToSegmentOld{T <: Number} <: Setup @deftype T
    D1
    H1
    D2
    H2
    σ
end
SegmentToSegment(old::SegmentToSegmentOld) = SegmentToSegment(D1=old.D1, H1=old.H1, D2=old.D2, H2=old.H2, σ=old.σ)

function precompute_z_weights(setup::SegmentToSegmentOld; params::Constants)
    @unpack D1, H1, D2, H2, σ = setup
    @unpack rb, line_points = params

    zt, wz = gausslegendre(sum(line_points)+1)   
    J1 = H1/2
    J2 = H2/2

    ζ1 = @. J1 * (zt + 1) + D1
    ζ2 = @. J2 * (zt + 1) + D2
    
    R̃ = [sqrt(σ^2 + (z2 - z1)^2) / rb for z1 in ζ1 for z2 in ζ2]
    w = [w1*w2 for w1 in wz for w2 in wz]
    return (R̃=R̃, w=w, J=J1*J2/H2)
end

function precompute_coefficients(setup::SegmentToSegmentOld; params::Constants, dp)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg = params

    R̃, wz, J = precompute_z_weights(setup, params=params)
    C = J * sqrt(m*π/2) / (2 * π^2 * rb * kg)

    P = zeros(n+1, n+1)
    R = zeros(n+1, length(R̃))
    M = zeros(n+1)

    for k in 0:n
        for s in 1:n+1
            P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    for (i, r̃) in enumerate(R̃)
        for k in 0:n
            R[k+1, i] = besselj(k+1/2, m * r̃) * imag((im)^k * exp(im*c*r̃)) / r̃^(3/2)
        end
    end

    M .= R * wz
    M .= P * M 
    return C .* M
end


function constant_integral(setup::SegmentToSegmentOld; params::Constants)
    @unpack kg, α = params
    @unpack D1, H1, D2, H2, σ = setup
    I(d) = sqrt(σ^2+d^2) + d * log(sqrt(σ^2+d^2) - d)
    1/(4π*kg*H2) * (I(D1+H1-D2-H2) + I(D1-D2) - I(D1-H2-D2) - I(D1+H1-D2))
end

function has_heatwave_arrived(setup::SegmentToSegmentOld; params::Constants, t)
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
