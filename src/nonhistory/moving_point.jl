@with_kw struct MovingPointToPoint{T <: Number} <: Setup @deftype T
    # Assuming movement of medium in the x-direction
    x   # x difference between source and evaluation point
    σ   # distance in the y, z plane between source and evaluation point
    v   # speed
    r = sqrt(x^2+σ^2)
end
gettype(s::MovingPointToPoint) = eltype(s.x)

function integral_formula(setup::MovingPointToPoint, params::Constants, fx, w, q, I_c) 
    @unpack x, v = setup
    @unpack α, Δt = params

    I = dot(fx, w)
    sig = sign(I)
    Ilog = log(abs(I)) + v/(4α) * (2x - v*Δt)
    #@show I
    #@show v/(4α) * (2x - v*Δt)
    #@show q * I_c
    sig * exp(Ilog) + q * I_c
end

function f_evolve_1!(setup::MovingPointToPoint, fx, x, Ct, q, params::Constants)
    @unpack rb, α, Δt = params
    @unpack v = setup
    ṽ = v*rb/α
    @. fx = Ct * (fx - q * x / (x^2 + ṽ^2/4))
    #@show maximum(abs.(Ct)), minimum(abs.(Ct))
    #@show maximum(abs.(fx)), minimum(abs.(fx))
end
function f_evolve_2!(setup::MovingPointToPoint, fx, x, q, params::Constants) 
    @unpack rb, α, Δt = params
    @unpack v = setup 
    ṽ = v*rb/α
    @. fx = exp(-v^2*Δt/(4α)) * fx + q * x / (x^2 + ṽ^2/4)
end
function f_guess(setup::MovingPointToPoint, params::Constants) 
    @unpack Δt, rb, α = params
    @unpack v = setup
    ṽ = v*rb / α
    #f(z) = exp(-10*z^2*Δt) * (1 - exp(-z^2*Δt)) * z / (z^2 + ṽ^2/4)
    f(z) = (1 - exp(-z^2*Δt)) * z / (z^2 + ṽ^2/4)
    f
end

function precompute_coefficients(setup::MovingPointToPoint; params::Constants, dp::DiscretizationParameters, containers::ComputationContainers)
    @unpack m, c, n, xt, w = dp
    @unpack r, σ, x, v = setup
    @unpack rb, kg, α, Δt = params
    r̃ = r/rb
    Ω = m * r̃
   
    C = m * exp(im * r̃ * c) / (2 * π^2 * r * kg)
    w = [ sum([imag(C * (im)^k) * (2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(xt[s],k) for k =0:n]) * w[s] for s=1:n+1]
    #@show maximum(abs.(w)), minimum(abs.(w))
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
