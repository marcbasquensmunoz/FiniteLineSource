
@with_kw struct SegmentToPoint{T <: Number} <: Setup @deftype T
    D
    H
    r
    z
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

    mul!(M, P, R)

    return C .* M * wz
end

function constant_integral(setup::SegmentToPoint; params::Constants)
    @unpack kg, α = params
    @unpack D, H, z, r = setup
    1 / (4π * kg) * log((z-D + sqrt(r^2 + (z-D)^2))/(z-D-H + sqrt(r^2 + (z-D-H)^2)))
end

function precompute_parameters(setup::SegmentToPoint; prealloc::Preallocation, params::Constants)
    @unpack segment_limits, segment_points = params
    dps = @views [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1], segment_limits[2:end], segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps)) 
    v  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params, P=prealloc.P[i], R=prealloc.R[i], M=prealloc.M[i]) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)

    return Precomputation(x=x, fx=fx, v=v, I_c=I_c)
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

function segment_to_point_test(setup::SegmentToPoint; params::Constants, t) 
    @unpack kg, α = params
    @unpack D, H, z, r = setup
    quadgk(ζ -> 1/(4π*kg*sqrt(r^2+(z-ζ)^2)) * erfc(sqrt(r^2+(z-ζ)^2)/sqrt(4α*t)), D, H+D, atol=10^-16, order=20)
end