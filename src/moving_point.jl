@with_kw struct MovingPointToPoint{T <: Number} <: Setup @deftype T
    # Assuming movement in the x-direction
    x   # initial x position
    σ   # distance to the axis of movement from the origin
    v   # speed
    r = sqrt(x^2 + σ^2)
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

function precompute_coefficients(setup::MovingPointToPoint; dp, params::Constants)
    @unpack m, c, n, xt, w = dp
    @unpack r, x, v = setup
    @unpack rb, kg, α = params
    r̃ = r/rb

    Ω = m * r̃
    C = m * exp(im * r̃ * c) * exp(x*v/(2α)) / (2 * π^2 * r * kg)
    w = [ sum([imag(C * (im)^k) *(2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(xt[s],k) for k =0:dp.n]) * w[s] for s=1:dp.n+1]

    return w
end

function constant_integral(setup::MovingPointToPoint; params::Constants)
    @unpack x, r, σ, v = setup
    @unpack kg, α = params
    exp(v*(x-r) / (2α)) / (4 * π * r * kg)
end

function has_heatwave_arrived(setup::MovingPointToPoint; params::Constants, t)
    true
end
