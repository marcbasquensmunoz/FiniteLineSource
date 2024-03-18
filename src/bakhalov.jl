abstract type Setup end

@with_kw struct Constants{T <: Number} @deftype T
    rb = .1
    α  = 10^-6
    kg = 3.
    Δt = 3600.

    segment_limits::Vector{T} = [0., 0.01, 0.1, 0.5, 1., 3., 5., 7., 10.] 
    segment_points::Vector{Int} = [10, 25, 15, 15, 20, 15, 15, 15]
    line_points::Vector{Int} = [30]
    line_limits::Vector{T} = [0., 1.]
end

@with_kw struct Preallocation{T <: Number}
    P::Vector{Matrix{T}}
    R::Vector{Matrix{T}}
    M::Vector{Vector{T}}
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

function precompute_parameters(setup::Setup; prealloc::Preallocation, params::Constants)
    @unpack segment_limits, segment_points = params
    dps = @views [discretization_parameters(a,b,n) for (a,b,n) in zip(segment_limits[1:end-1], segment_limits[2:end], segment_points)] 
    
    x  = reduce(vcat, (dp.x for dp in dps)) 
    v  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params, P=prealloc.P[i], R=prealloc.R[i], M=prealloc.M[i]) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)

    return Precomputation(x=x, fx=fx, v=v, I_c=I_c)
end

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end
