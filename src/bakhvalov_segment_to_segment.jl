using Parameters
using FastGaussQuadrature
using SpecialFunctions
using LegendrePolynomials
using Bessels: besselj
using LinearAlgebra
using Cubature
using QuadGK

const threshold = 8. 

@with_kw struct SimulationParametersSegToSeg{T <: Real} @deftype T
    D1
    H1
    D2
    H2
    Δr
    rb = 0.1
    α = 10^-6
    kg = 3.
    Δt = 3600*24*30.
    segment_limits::Vector{T} = [0., 10.]
    segment_points::Vector{Int} = [100]
    line_points::Int
end

@with_kw struct Preallocations{T <: Real} @deftype T
    P::Vector{Matrix{T}}
    R::Vector{Matrix{T}}
    M::Vector{Matrix{T}}
end
function Preallocations(segment_points, line_points) 
    P = [zeros(i+1, i+1) for i in segment_points]
    R = [zeros(i+1, (line_points+1)^2) for i in segment_points]
    M = [zeros(i+1, (line_points+1)^2) for i in segment_points]
    Preallocations(P=P, R=R, M=M)
end

@with_kw struct Precomputation{T <: Real} @deftype T
    x::Vector{T}
    fx::Vector{T}
    v::Vector{T}
end

function has_heatwave_arrived(params::SimulationParametersSegToSeg; t) 
    @unpack D1, H1, D2, H2, Δr, α = params
    true
    #=if z >= D && z <= D + H
        r^2 / (4α*t) < threshold^2
    else 
        r^2 + min((D-z)^2, (D+H-z)^2) / (4α*t) < threshold^2
    end=#
end

function segment_to_segment(;params::SimulationParametersSegToSeg, t) 
    @unpack D1, H1, D2, H2, Δr, kg, α = params
    quadgk(z -> quadgk(ζ -> q[1]/H2/(4π*kg*sqrt(Δr^2+(z-ζ)^2)) * erfc(sqrt(Δr^2+(z-ζ)^2)/sqrt(4α*t)), D1, H1+D1, atol=10^-16, order=20)[1], D2, D2+H2, atol=10^-16, order=20)
end

function compute_integral_slow(q; params::SimulationParametersSegToSeg) 
    @unpack D1, H1, D2, H2, Δr, kg = params
    hcubature(z -> q / H2 / (4π * kg) * log((z[1]-D2 + sqrt(Δr^2 + (z[1]-D2)^2))/(z[1]-D2-H2 + sqrt(Δr^2 + (z[1]-D2-H2)^2))), D1, D1+H1)[1]
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

function precompute_z_weights(params::SimulationParametersSegToSeg)
    @unpack D1, H1, D2, H2, Δr, rb, α, line_points = params

    zt, wz = gausslegendre(line_points+1)   
    J1 = H1/2
    J2 = H2/2

    ζ1 = @. J1 * (zt + 1) + D1
    ζ2 = @. J2 * (zt + 1) + D2
    
    R̃ = [sqrt(Δr^2 + (z2 - z1)^2) / rb for z1 in ζ1 for z2 in ζ2]
    w = [w1*w2 for w1 in wz for w2 in wz]

    return (R̃=R̃, w=w, J=J1*J2)
end

function precompute_coefficients(params::SimulationParametersSegToSeg; dp, P, R, M)
    @unpack m, c, n, xt, w = dp
    @unpack D1, H1, D2, H2, Δr, rb, kg, α = params

    R̃, wz, J = precompute_z_weights(params)
    C = J / H2 * sqrt(m*π/2) / (2 * π^2 * rb * kg)

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

function precompute_parameters(params::SimulationParametersSegToSeg; prealloc::Preallocations)
    @unpack segment_limits, segment_points = params
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1],segment_limits[2:end],segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps))
    v  = reduce(vcat, [precompute_coefficients(params, dp=dp, P=prealloc.P[i], R=prealloc.R[i], M=prealloc.M[i]) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))

    return Precomputation(x=x, fx=fx, v=v)
end

function compute_integral_throught_history!(I, q; precomp::Precomputation, params::SimulationParametersSegToSeg)
    @unpack D1, H1, D2, H2, Δr, rb, kg, α, Δt = params
    @unpack x, fx, v = precomp
 
    Δt̃ = α*Δt/rb^2
    Ct = @. exp(-x^2*Δt̃)
    for k in eachindex(I)
        @views @inbounds f_evolve_1!(fx, x, Ct, q[k])
        if has_heatwave_arrived(params, t=Δt*k)
            @views @inbounds I[k] = compute_integral_oscillatory(fx, v) + compute_integral_slow(q[k], params=params)
        else 
            I[k] = 0.
        end
        @views @inbounds f_evolve_2!(fx, x, q[k])    
    end
    
    return nothing
end

### Mean effect of bh 1 on bh 2

## Example
segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 5., 7., 10.] 
segment_points = 2 .* [10, 25, 15, 15, 20, 15, 15, 15]
params = SimulationParametersSegToSeg(D1=0., H1=10., D2=0., H2=10., Δr=1., segment_limits=segment_limits, segment_points=segment_points, line_points=10)

q = ones(1)
nt = length(q)
I = zeros(nt)

prealloc = Preallocations(segment_points, params.line_points) 

@time precomp = precompute_parameters(params, prealloc=prealloc)
@time compute_integral_throught_history!(I, q, precomp=precomp, params=params)
I                                            
res = segment_to_segment(params=params, t=params.Δt) 
err = res[1] - I[1]
