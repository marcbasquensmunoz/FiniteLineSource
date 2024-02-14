@with_kw struct SegmentToSegment{T <: Number} <: Setup @deftype T
    D1
    H1
    D2
    H2
    r
end

function precompute_z_weights(setup::SegmentToSegment; params::Constants)
    @unpack D1, H1, D2, H2, r = setup
    @unpack rb, line_points = params

    zt, wz = gausslegendre(line_points+1)   
    J1 = H1/2
    J2 = H2/2

    ζ1 = @. J1 * (zt + 1) + D1
    ζ2 = @. J2 * (zt + 1) + D2
    
    R̃ = [sqrt(r^2 + (z2 - z1)^2) / rb for z1 in ζ1 for z2 in ζ2]
    w = [w1*w2 for w1 in wz for w2 in wz]

    return (R̃=R̃, w=w, J=J1*J2/H2)
end

function precompute_coefficients(setup::SegmentToSegment; params::Constants, dp, P, R, M)
    @unpack m, c, n, xt, w = dp
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

function constant_integral(setup::SegmentToSegment; params::Constants)
    @unpack kg, α = params
    @unpack D1, H1, D2, H2, r = setup
    hcubature(z -> 1 / (4π * kg * H2) * log((z[1]-D2 + sqrt(r^2 + (z[1]-D2)^2))/(z[1]-D2-H2 + sqrt(r^2 + (z[1]-D2-H2)^2))), D1, D1+H1)[1]
end

function precompute_parameters(setup::SegmentToSegment; prealloc::Preallocation, params::Constants)
    @unpack segment_limits, segment_points = params
    dps = @views [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1], segment_limits[2:end], segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps)) 
    v  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params, P=prealloc.P[i], R=prealloc.R[i], M=prealloc.M[i]) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)

    return Precomputation(x=x, fx=fx, v=v, I_c=I_c)
end


function segment_to_segment_test(setup::SegmentToSegment; params::Constants, t) 
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

