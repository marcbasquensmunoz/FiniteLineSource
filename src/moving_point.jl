@with_kw struct MovingPointToPoint{T <: Number} <: Setup @deftype T
    # Assuming movement in the x-direction
    x   # initial x position
    σ   # distance to the axis of movement from the origin
    v   # speed
    r = sqrt(x^2 + σ^2)
end

function Preallocation(::MovingPointToPoint, params::Constants) 
    @unpack segment_points = params
    Preallocation(P=[zeros(0, 0) for i in segment_points], R=[zeros(0, 0) for i in segment_points], M=[zeros(0) for i in segment_points])
end

function f_evolve_1!(setup::MovingPointToPoint, fx, x, Ct, q, params::Constants)
    @unpack rb, α, Δt = params
    @unpack v = setup

    @. fx = exp(-v^2*Δt/(4α)) * Ct * (fx - q * x / (x^2 + (v * rb / α)^2/4))
end
function f_evolve_2!(setup::MovingPointToPoint, fx, x, q, params::Constants) 
    @unpack rb, α = params
    ṽ = setup.v*rb / α
    @. fx = fx + q * x / (x^2 + ṽ^2/4)
end

function precompute_coefficients(setup::MovingPointToPoint; dp, params::Constants, P, R, M)
    @unpack m, c, n, xt, w = dp
    @unpack r, σ, v = setup
    @unpack rb, kg, α = params
    r̃ = r/rb

    Ω = dp.m * r̃
    C = dp.m * exp(im * r̃ * dp.c) * exp(r*v/(2α)) / (2 * π^2 * r * kg)
    w = [ sum([imag(C * (im)^k) *(2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(dp.xt[s],k) for k =0:dp.n]) * dp.w[s] for s=1:dp.n+1]

    return w
end

function constant_integral(setup::MovingPointToPoint; params::Constants)
    @unpack r, σ, v = setup
    @unpack kg = params
    1 / (4 * π * r * kg)
end

function has_heatwave_arrived(setup::MovingPointToPoint; params::Constants, t)
    @unpack r, v = setup
    @unpack α = params
    threshold = 10^-20

    (r+t*v)/sqrt(4*t*α) > erfcinv(threshold * exp(-r*abs(v)/α))
end

function analytical_test(setup::MovingPointToPoint; params::Constants, t)
    @unpack r, σ, v = setup
    @unpack α, kg = params
    moving_point_step_response(t, r, σ, v, α, kg)
end
