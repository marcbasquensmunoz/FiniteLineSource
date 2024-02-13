using Parameters
using FastGaussQuadrature
using SpecialFunctions
using LegendrePolynomials
using Bessels: besselj
using LinearAlgebra
using QuadGK

const threshold = 8. 

@with_kw struct SimulationParameters{T <: Real} @deftype T
    D
    H
    rb = 0.1
    α = 10^-6
    kg = 3.
    Δt = 3600*24*30.
    segment_limits::Vector{T} = [0., 10.]
    segment_points::Vector{Int} = [100]
    line_points::Int
end

@with_kw struct Evaluation{T <: Real} @deftype T
    z
    r
end

@with_kw struct Preallocations{T <: Real} @deftype T
    P::Vector{Matrix{T}}
    R::Vector{Matrix{T}}
    M::Vector{Matrix{T}}
end
function Preallocations(segment_points, line_points) 
    P = [zeros(i+1, i+1) for i in segment_points]
    R = [zeros(i+1, line_points+1) for i in segment_points]
    M = [zeros(i+1, line_points+1) for i in segment_points]
    Preallocations(P=P, R=R, M=M)
end

@with_kw struct Precomputation{T <: Real} @deftype T
    x::Vector{T}
    fx::Vector{T}
    v::Vector{T}
end

function has_heatwave_arrived(ev::Evaluation, params::SimulationParameters; t) 
    @unpack D, H, α = params
    @unpack z, r = ev
    if z >= D && z <= D + H
        r^2 / (4α*t) < threshold^2
    else 
        r^2 + min((D-z)^2, (D+H-z)^2) / (4α*t) < threshold^2
    end
end

function segment_to_point(;ev::Evaluation, params::SimulationParameters, t) 
    @unpack D, H, rb, kg, α = params
    @unpack z, r = ev
    quadgk(ζ -> q[1]/(4π*kg*sqrt(r^2+(z-ζ)^2)) * erfc(sqrt(r^2+(z-ζ)^2)/sqrt(4α*t)), D, H+D, atol=10^-16, order=20)
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

function precompute_z_weights(;params::SimulationParameters, ev::Evaluation)
    @unpack D, H, rb, α, line_points = params
    @unpack r, z = ev

    z1 = D
    z2 = D+H
    zt, wz = gausslegendre(line_points+1)   
    J = (z2-z1)/2
    for (i, x) in enumerate(zt)
        zt[i] = @. sqrt(r^2 + (z - (J * (x + 1) + z1))^2) / rb
    end
    return (R=zt, wz=wz, J=J)
end

function precompute_coefficients(dp; params::SimulationParameters, ev::Evaluation, P, R, M)
    @unpack m, c, n, xt, w = dp
    @unpack D, H, rb, kg, α = params
    @unpack r, z = ev

    R̃, wz, J = precompute_z_weights(params=params, ev=ev)
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

function precompute_parameters(;ev::Evaluation, params::SimulationParameters, prealloc::Preallocations)
    @unpack segment_limits, segment_points = params
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1],segment_limits[2:end],segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps))
    v  = reduce(vcat, [precompute_coefficients(dp, params=params, ev=ev, P=prealloc.P[i], R=prealloc.R[i], M=prealloc.M[i]) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))

    return Precomputation(x=x, fx=fx, v=v)
end

function compute_integral_throught_history!(I, q; precomp::Precomputation, ev::Evaluation, params::SimulationParameters)
    @unpack D, H, rb, kg, α, Δt = params
    @unpack x, fx, v = precomp
 
    Δt̃ = α*Δt/rb^2
    Ct = @. exp(-x^2*Δt̃)
    for k in eachindex(I)
        @views @inbounds f_evolve_1!(fx, x, Ct, q[k])
        if has_heatwave_arrived(ev, params, t=Δt*k)
            @views @inbounds I[k] = compute_integral_oscillatory(fx, v) + compute_integral_slow(q[k], params=params, ev=ev)
        else 
            I[k] = 0.
        end
        @views @inbounds f_evolve_2!(fx, x, q[k])    
    end
    
    return nothing
end


## Example
segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 5., 7., 10.] 
segment_points = 2 .* [10, 25, 15, 15, 20, 15, 15, 15]
params = SimulationParameters(D=0., H=50., segment_limits=segment_limits, segment_points=segment_points, line_points=50)

q = ones(1)
nt = length(q)
I = zeros(nt)
ev = Evaluation(z=25., r=1.)

prealloc = Preallocations(segment_points, params.line_points) 

@time precomp = precompute_parameters(ev=ev, params=params, prealloc=prealloc)
@time compute_integral_throught_history!(I, q, precomp=precomp, params=params, ev=ev)
I                                            
res = segment_to_point(ev=ev, params=params, t=params.Δt) 
err = res[1] - I[1]
