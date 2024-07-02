@with_kw struct MovingPointToPoint{T <: Number} <: Setup @deftype T
    # Assuming movement of medium in the x-direction
    x   # x difference between source and evaluation point
    σ   # distance in the y, z plane between source and evaluation point
    v   # speed
    r = sqrt(x^2+σ^2)
end
gettype(s::MovingPointToPoint) = eltype(s.x)

function f_evolve_1!(setup::MovingPointToPoint, fx, x, Ct, q, params::Constants)
    @unpack rb, α, Δt = params
    @unpack v = setup
    ṽ = v*rb/α
    @. fx = exp(-v^2*Δt/(4α)) * Ct * (fx - q * x / (x^2 + ṽ^2/4))
end
function f_evolve_2!(setup::MovingPointToPoint, fx, x, q, params::Constants) 
    @unpack rb, α = params
    @unpack v = setup 
    ṽ = v*rb/α
    @. fx = fx + q * x / (x^2 + ṽ^2/4)
end
function f_guess(setup::MovingPointToPoint, params::Constants) 
    @unpack Δt, rb, α = params
    @unpack v = setup
    ṽ = v*rb / α
    f(z) = exp(-10*z^2*Δt) * (1 - exp(-z^2*Δt)) * z / (z^2 + ṽ^2/4)
    f
end

function precompute_coefficients(setup::MovingPointToPoint; dp, params::Constants)
    @unpack m, c, n, xt, w = dp
    @unpack r, σ, x, v = setup
    @unpack rb, kg, α = params
    r̃ = r/rb
    Ω = m * r̃
   
    E = exp(x*v/(2α)) #x*v/(2α) < 710 ? exp(x*v/(2α)) : 1 # Avoid NaN
    C = m * exp(im * r̃ * c) * E / (2 * π^2 * r * kg)
    w = [ sum([imag(C * (im)^k) * (2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(xt[s],k) for k =0:n]) * w[s] for s=1:n+1]

    return w
end

function constant_integral(setup::MovingPointToPoint; params::Constants)
    @unpack x, r, v = setup
    @unpack kg, α = params
    exp(-(r-x)*v / (2α)) / (4 * π * r * kg)
end

function has_heatwave_arrived(setup::MovingPointToPoint; params::Constants, t)
    @unpack x, r, v = setup
    @unpack α = params
    exp(v*(x-r)/(2α)) * erfc((r-t*v)/sqrt(4α*t)) > 10^-20 || exp(v*(x+r)/(2α)) * erfc((r+t*v)/sqrt(4α*t)) > 10^-20
end
