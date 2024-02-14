
@with_kw struct PointToPoint{T <: Number} <: Setup @deftype T
    r
end

function precompute_coefficients(setup::PointToPoint; dp, params::Constants)
    @unpack m, c, n, xt, w = dp
    @unpack r = setup
    @unpack rb = params

    Ω = dp.m * r / rb
    C = dp.m * exp(im * r / rb * dp.c) / (2 * π^2 * r * kg)
    v = [ sum([imag(C * (im)^k) *(2k+1) * sqrt(π/(2Ω)) * besselj(k+1/2, Ω)*Pl.(dp.xt[s],k) for k =0:dp.n]) * dp.w[s] for s=1:dp.n+1]

    return v
end

function constant_integral(setup::PointToPoint; params::Constants)
    @unpack r = setup
    @unpack kg = params
    π/2 / (2 * π^2 * r * kg)
end

function precompute_parameters(setup::PointToPoint; prealloc::Preallocation, params::Constants)
    @unpack segment_limits, segment_points = params
    dps = @views [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1], segment_limits[2:end], segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps)) 
    v  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)

    return Precomputation(x=x, fx=fx, v=v, I_c=I_c)
end

function has_heatwave_arrived(setup::PointToPoint; params::Constants, t)
    @unpack r = setup
    @unpack α = params
    threshold = 8
    r^2 / (2α*t) < threshold^2
end

function point_to_point_test(setup::PointToPoint; params::Constants, t)
    @unpack r = setup
    @unpack α, kg = params
    erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)
end
