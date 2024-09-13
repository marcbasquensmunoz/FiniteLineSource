@with_kw struct MovingSegmentToPoint{T <: Number} <: Setup @deftype T
    # Assuming movement in the x-direction
    x   
    y
    z
    v
    D
    H
    σ = sqrt(x^2 + y^2)
end

function f_evolve_1!(setup::MovingSegmentToPoint, fx, x, Ct, q, params::Constants)
    @unpack rb, α, Δt = params
    @unpack v = setup

    @. fx = exp(-v^2*Δt/(4α)) * Ct * (fx - q * x / (x^2 + (v * rb / α)^2/4))
end
function f_evolve_2!(setup::MovingSegmentToPoint, fx, x, q, params::Constants) 
    @unpack rb, α = params
    ṽ = setup.v*rb / α
    @. fx = fx + q * x / (x^2 + ṽ^2/4)
end

function precompute_z_weights(setup::MovingSegmentToPoint; params::Constants)
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

function precompute_coefficients(setup::MovingSegmentToPoint; params::Constants, dp::DiscretizationParameters, containers::ComputationContainers)
    @unpack m, c, n, xt, w = dp
    @unpack rb, kg, α = params
    @unpack x, v = setup

    R̃, wz = precompute_z_weights(setup, params=params)
    C = m * exp(x*v/(2α)) / (2 * π^2 * kg)

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

function constant_integral(setup::MovingSegmentToPoint; params::Constants)
    @unpack x, z, σ, v, D, H = setup
    @unpack kg, α = params
    quadgk(ζ -> exp(v*(x-sqrt(σ^2 + (z-ζ)^2)) / (2α)) / (4 * π * sqrt(σ^2 + (z-ζ)^2) * kg), D, D+H)[1]
end

function has_heatwave_arrived(setup::MovingSegmentToPoint; params::Constants, t)
    true
end
