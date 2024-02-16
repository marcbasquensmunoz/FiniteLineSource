
@with_kw struct SegmentToPoint{T <: Number} <: Setup @deftype T
    D
    H
    r
    z
end

function Preallocation(::SegmentToPoint, params::Constants) 
    @unpack segment_points, line_points = params
    P = [zeros(i+1, i+1) for i in segment_points]
    R = [zeros(i+1, line_points+1) for i in segment_points]
    M = [zeros(i+1) for i in segment_points]
    Preallocation(P=P, R=R, M=M)
end

function precompute_z_weights(setup::SegmentToPoint; params::Constants)
    @unpack rb, α, line_points = params
    @unpack D, H, r, z = setup

    zt, wz = gausslegendre(line_points+1)   
    J = H/2
    for (i, x) in enumerate(zt)
        zt[i] = sqrt(r^2 + (z - (J * (x + 1) + D))^2) / rb
    end
    return (R=zt, wz=wz, J=J)
end

function precompute_coefficients(setup::SegmentToPoint; params::Constants, dp, P, R, M)
    @unpack m, c, n, xt, w = dp
    @unpack D, H, r, z = setup
    @unpack rb, kg = params

    R̃, wz, J = precompute_z_weights(setup, params=params)
    C = J * sqrt(m*π/2) / (2 * π^2 * rb * kg)

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

function constant_integral(setup::SegmentToPoint; params::Constants)
    @unpack kg, α = params
    @unpack D, H, z, r = setup
    1 / (4π * kg) * log((z-D + sqrt(r^2 + (z-D)^2))/(z-D-H + sqrt(r^2 + (z-D-H)^2)))
end

function has_heatwave_arrived(setup::SegmentToPoint; params::Constants, t)
    @unpack α = params
    @unpack D, H, z, r = setup
    threshold = 8

    if z >= D && z <= D + H
        r^2 / (4α*t) < threshold^2
    else 
        r^2 + min((D-z)^2, (D+H-z)^2) / (4α*t) < threshold^2
    end
end

function analytical_test(setup::SegmentToPoint; params::Constants, t) 
    @unpack kg, α = params
    @unpack D, H, z, r = setup
    quadgk(ζ -> 1/(4π*kg*sqrt(r^2+(z-ζ)^2)) * erfc(sqrt(r^2+(z-ζ)^2)/sqrt(4α*t)), D, H+D, atol=10^-16, order=20)
end