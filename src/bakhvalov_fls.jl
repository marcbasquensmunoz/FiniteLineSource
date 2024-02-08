using Parameters
using FastGaussQuadrature
using SpecialFunctions
using LegendrePolynomials
using Bessels
using LinearAlgebra
using QuadGK

@with_kw struct SimulationParameters{T <: Real} @deftype T
    D
    H
    rb = 0.1
    α = 10^-6
    kg = 3.
    segment_limits::Vector{T} = [0., 10.]
    segment_points::Vector{Int} = [100]
end

@with_kw struct Evaluation{T <: Real} @deftype T
    Δt
    z
    r
end

@with_kw struct Precomputation{T <: Real} @deftype T
    x::Vector{T}
    fx::Vector{T}
    v::Vector{T}
end

function segment_to_point(;ev::Evaluation, params::SimulationParameters) 
    @unpack D, H, rb, kg, α = params
    @unpack Δt, z, r = ev
    quadgk(ζ -> q[1]/(4π*kg*sqrt(r^2+(z-ζ)^2)) * erfc(sqrt(r^2+(z-ζ)^2)/sqrt(4α*Δt)), D, H+D)
end

function compute_integral_slow(q; params::SimulationParameters, ev::Evaluation) 
    @unpack D, H, kg = params
    @unpack z, r = ev
    q / (4π * kg) * log((z-D + sqrt(r^2 + (z-D)^2))/(z-D-H + sqrt(r^2 + (z-D-H)^2)))
end
compute_integral_oscillatory(fx, v) = dot(fx, v)

function f_evolve_1!(fx, x, Ct, q)
    @. fx = Ct * (fx - q/x)
end

function f_evolve_2!(fx, x, q)
    @. fx = fx + q / x
end

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end

function precompute_z_weights(Ns; params::SimulationParameters, r, z)
    @unpack D, H, rb = params
    zt, wz = gausslegendre(Ns+1)   
    ζ = @. H/2 * (zt + 1) + D
    R = @. sqrt(r^2 + (z - ζ)^2) / rb
    return (R=R, wz=wz)
end

function precompute_coefficients(dp; params::SimulationParameters, ev::Evaluation, Ns=1000)
    @unpack m, c, n, xt, w = dp
    @unpack D, H, rb, kg = params
    @unpack r, z = ev
    C = H/2 * sqrt(m*π/2) / (2 * π^2 * rb * kg)
    R, wz = precompute_z_weights(Ns, params=params, r=r, z=z)

    P = reduce(hcat, [[w[s] * (2k+1) * Pl.(xt[s], k) for s in 1:n+1] for k in 0:n])
    R = reduce(hcat, [[Bessels.besselj(k+1/2, m * rr) * imag((im)^k * exp(im*c*rr)) / rr^(3/2) for k in 0:n] for rr in R])

    return C .* P * R * wz
end

function precompute_parameters(;ev::Evaluation, params::SimulationParameters)
    @unpack segment_limits, segment_points = params
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1],segment_limits[2:end],segment_points)] 
    fps = [precompute_coefficients(dp, params=params, ev=ev) for dp in dps]
    fxs = [zeros(dp.n+1) for dp in dps]
    
    x  = reduce(vcat, (dp.x for dp in dps))
    fx = reduce(vcat, fxs)
    v  = reduce(vcat, (fp for fp in fps))

    return Precomputation(x=x, fx=fx, v=v)
end

function compute_integral_throught_history!(I, q; precomp::Precomputation, ev::Evaluation, params::SimulationParameters)
    @unpack D, H, rb, kg, α = params
    @unpack Δt = ev
    @unpack x, fx, v = precomp
    Δt̃ = α*Δt/rb^2
    Ct = @. exp(-x^2*Δt̃)
    for k in eachindex(I)
        @views @inbounds f_evolve_1!(fx, x, Ct, q[k])
        @views @inbounds I[k] = compute_integral_oscillatory(fx, v) + compute_integral_slow(q[k], params=params, ev=ev)
        @views @inbounds f_evolve_2!(fx, x, q[k])    
    end
    
    return nothing
end


## Example
segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 5., 7., 10.] 
segment_points = [10, 25, 15, 15, 20, 15, 15, 15]
params = SimulationParameters(D=0., H=20., segment_limits=segment_limits, segment_points=segment_points)
ev = Evaluation(Δt = 3600*24*365*5., z = 10., r = 0.5)

q = [1.]
I = zeros(length(q))

precomp = precompute_parameters(ev=ev, params=params)
compute_integral_throught_history!(I, q, precomp=precomp, params=params, ev=ev)
I                                      # Compare this result 

segment_to_point(ev=ev, params=params) # Compare to this
