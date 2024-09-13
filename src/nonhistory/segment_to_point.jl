
@with_kw struct SegmentToPoint{T <: Number} <: Setup @deftype T
    D
    H
    σ
    z
end

function precompute_z_weights(setup::SegmentToPoint; params::Constants)
    @unpack rb, α, line_points, line_limits = params
    @unpack D, H, σ, z = setup

    R = zeros(0)
    wz = zeros(0)

    limits = line_limits .* H .+ D

    for (a, b, n) in zip(limits[1:end-1], limits[2:end], line_points) 
        @unpack x, m, w = discretization_parameters(a, b, n)
        append!(R, @. sqrt(σ^2 + (z - x)^2) / rb)
        append!(wz, m .* w)
    end

    return (R=R, wz=wz)
end

function precompute_coefficients(setup::SegmentToPoint; params::Constants, dp, containers::ComputationContainers)
    @unpack m, c, n, xt, w = dp
    @unpack D, H, z = setup
    @unpack rb, kg = params

    R̃, wz = precompute_z_weights(setup, params=params)
    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

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

function constant_integral(setup::SegmentToPoint; params::Constants)
    @unpack kg, α = params
    @unpack D, H, z, σ = setup
    1 / (4π * kg) * log((z-D + sqrt(σ^2 + (z-D)^2))/(z-D-H + sqrt(σ^2 + (z-D-H)^2)))
end

function has_heatwave_arrived(setup::SegmentToPoint; params::Constants, t)
    @unpack α = params
    @unpack D, H, z, σ = setup
    threshold = 8

    if z >= D && z <= D + H
        σ^2 / (4α*t) < threshold^2
    else 
        σ^2 + min((D-z)^2, (D+H-z)^2) / (4α*t) < threshold^2
    end
end
