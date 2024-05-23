abstract type Setup end
gettype(::Setup) = Float64

@with_kw struct Constants{T <: Number} @deftype T
    rb = .1
    α  = 10^-6
    kg = 3.
    Δt = 3600.

    b = 10.
    line_points::Vector{Int} = [30, 30, 30]
    line_limits::Vector{T} = [0., 0.4, 0.6, 1.]
end

@with_kw struct Precomputation{T <: Number}
    x::Vector{T}
    fx::Vector{T}
    w::Vector{T}
    I_c::T
end

integral_formula(::Setup, ::Constants, fx, w, q, I_c) = dot(fx, w) + q * I_c
f_evolve_1!(::Setup, fx, x, Ct, q, params::Constants) = @. fx = Ct * (fx - q/x)
f_evolve_2!(::Setup, fx, x, q, params::Constants)     = @. fx = fx + q / x
function f_guess(setup::Setup, params::Constants) 
    @unpack Δt = params
    f(z) = exp(-10*z^2*Δt) * (1 - exp(-z^2*Δt)) / z
    f
end

function compute_integral_throught_history!(setup::Setup; I, q, precomp::Precomputation, params::Constants)
    @unpack rb, α, Δt = params
    @unpack x, fx, w, I_c = precomp

    Δt̃ = α*Δt/rb^2
    Ct = @. exp(-x^2*Δt̃)
    for k in eachindex(I)
        @inbounds f_evolve_1!(setup, fx, x, Ct, q[k], params)
        if has_heatwave_arrived(setup, params=params, t=Δt*k)
            @inbounds I[k] = integral_formula(setup, params, fx, w, q[k], I_c)
        else 
            @inbounds I[k] = 0.
        end
        @inbounds f_evolve_2!(setup, fx, x, q[k], params)    
    end
    
    return nothing
end

function precompute_parameters(setup::Setup; params::Constants, n = 20, n_tot = [0])
    @unpack Δt, b = params
    type = gettype(setup)

    segments = adaptive_gk_segments(f_guess(setup, params), convert(type, 0.), convert(type, b))
    dps = @views [discretization_parameters(s.a, s.b, n) for s in segments]
    x  = reduce(vcat, (dp.x for dp in dps)) 
    w  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params) for (i, dp) in enumerate(dps)])
    fx = zeros(type, sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)
    n_tot[1] = length(x)

    return Precomputation(x=x, fx=fx, w=w, I_c=I_c)
end

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end
