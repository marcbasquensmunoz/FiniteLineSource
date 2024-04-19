
@with_kw struct PointToPoint{T <: Number} <: Setup @deftype T
    r
end

function Preallocation(::PointToPoint, params::Constants) 
    @unpack segment_points = params
    Preallocation(P=[zeros(0, 0) for i in segment_points], R=[zeros(0, 0) for i in segment_points], M=[zeros(0) for i in segment_points])
end

function precompute_coefficients(setup::PointToPoint; dp, params::Constants)
    @unpack m, c, n, xt, w = dp
    @unpack r = setup
    @unpack rb, kg = params

    Ω = m * r / rb
    C = m * exp(im * r / rb * c) / (2 * π^2 * r * kg)
    v = [ sum([imag(C * (im)^k) * (2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(xt[s],k) for k =0:n]) * w[s] for s=1:n+1]

    return v
end

function constant_integral(setup::PointToPoint; params::Constants)
    @unpack r = setup
    @unpack kg = params
    π/2 / (2 * π^2 * r * kg)
end

function has_heatwave_arrived(setup::PointToPoint; params::Constants, t)
    @unpack r = setup
    @unpack α = params
    threshold = 8
    r^2 / (2α*t) < threshold^2
end

function analytical_test(setup::PointToPoint; params::Constants, t)
    @unpack r = setup
    @unpack α, kg = params
    erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)
end
