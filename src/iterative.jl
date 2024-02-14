abstract type Setup end

@with_kw struct Constants{T <: Number} @deftype T
    rb = .1
    α  = 10^-6
    kg = 3.
    Δt = 3600.

    segment_limits::Vector{T} = [0., 0.01, 0.1, 0.5, 1., 3., 5., 7., 10.] 
    segment_points::Vector{Int} = 2 .* [10, 25, 15, 15, 20, 15, 15, 15]
    line_points::Int = 30
end

@with_kw struct Preallocation{T <: Number}
    P::Vector{Matrix{T}}
    R::Vector{Matrix{T}}
    M::Vector{Matrix{T}}
end
function Preallocation(::SegmentToPoint, params::Constants) 
    @unpack segment_points, line_points = params
    P = [zeros(i+1, i+1) for i in segment_points]
    R = [zeros(i+1, line_points+1) for i in segment_points]
    M = [zeros(i+1, line_points+1) for i in segment_points]
    Preallocation(P=P, R=R, M=M)
end
function Preallocation(::SegmentToSegment, params::Constants) 
    @unpack segment_points, line_points = params
    P = [zeros(i+1, i+1) for i in segment_points]
    R = [zeros(i+1, (line_points+1)^2) for i in segment_points]
    M = [zeros(i+1, (line_points+1)^2) for i in segment_points]
    Preallocation(P=P, R=R, M=M)
end
function Preallocation(::PointToPoint, params::Constants) 
    Preallocation(P=[zeros(0, 0)], R=[zeros(0, 0)], M=[zeros(0, 0)])
end

@with_kw struct Precomputation{T <: Number}
    x::Vector{T}
    fx::Vector{T}
    v::Vector{T}
    I_c::T
end

exponential_integral(fx, v) = dot(fx, v)
f_evolve_1!(fx, x, Ct, q) = @. fx = Ct * (fx - q/x)
f_evolve_2!(fx, x, q)     = @. fx = fx + q / x

function compute_integral_throught_history!(setup::Setup; I, q, precomp::Precomputation, params::Constants)
    @unpack rb, α, Δt = params
    @unpack x, fx, v, I_c = precomp
 
    Δt̃ = α*Δt/rb^2
    Ct = @. exp(-x^2*Δt̃)
    for k in eachindex(I)
        @inbounds f_evolve_1!(fx, x, Ct, q[k])
        if has_heatwave_arrived(setup, params=params, t=Δt*k)
            @inbounds I[k] = exponential_integral(fx, v) + q[k] * I_c
        else 
            @inbounds I[k] = 0.
        end
        @inbounds f_evolve_2!(fx, x, q[k])    
    end
    
    return nothing
end

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end



q = ones(1)
nt = length(q)
I = zeros(nt)
setup = PointToPoint(r=5.)
params = Constants(Δt = 3600*24*30.)
prealloc = Preallocation(setup, params) 

@time precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
@time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
I                                            

res = point_to_point_test(setup, params=params, t=params.Δt) 
err = res[1] - I[1]