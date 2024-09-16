
@with_kw struct Constants{T <: Number} @deftype T
    rb = .1
    α  = 10^-6
    kg = 3.
    Δt = 3600.

    b = 10.
    line_points::Vector{Int} = [30, 30, 30]
    line_limits::Vector{T} = [0., 0.4, 0.6, 1.]
end

@with_kw mutable struct Precomputation{T <: Number}
    x::Vector{T}
    fx::Vector{T}
    w::Vector{T}
    I_c::T
end

@with_kw struct DiscretizationParameters{T <: Number} @deftype T
    xt::Vector{T}
    w::Vector{T}
    x::Vector{T}
    m
    c
    n::Int
end

function DiscretizationParameters(a, b, n::Int) 
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    DiscretizationParameters(xt=xt, w=w, x=x, m=m, c=c, n=n)
end

struct EmptyContainer <: ComputationContainers end

initialize_containers(::Setup, dps) = ([EmptyContainer()], ones(Int64, length(dps)))
