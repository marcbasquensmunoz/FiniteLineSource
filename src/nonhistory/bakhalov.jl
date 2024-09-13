

integral_formula(::Setup, ::Constants, fx, w, q, I_c) = dot(fx, w) + q * I_c
f_evolve_1!(::Setup, fx, x, Ct, q, params::Constants) = @. fx = Ct * (fx - q/x)
f_evolve_2!(::Setup, fx, x, q, params::Constants)     = @. fx = fx + q / x
function f_guess(::Setup, params::Constants) 
    @unpack Δt, α, rb = params
    Δt̃ = α*Δt/rb^2
    T = [1, 24*7, 24*365, 24*365*20]
    w = 10 .^ collect(0:length(T)-1)
    f(z) = sum(w .* exp.(-T*z^2*Δt̃)) * (1 - exp(-z^2*Δt̃)) / z
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
    dps = @views [DiscretizationParameters(s.a, s.b, n) for s in segments]
    containers, map = initialize_containers(setup, dps)    
    x  = reduce(vcat, (dp.x for dp in dps))
    w  = reduce(vcat, [precompute_coefficients(setup, dp=dp, params=params, containers=containers[map[i]]) for (i, dp) in enumerate(dps)])
    fx = zeros(type, sum([dp.n+1 for dp in dps]))
    I_c = constant_integral(setup, params=params)
    n_tot[1] = length(x)

    return Precomputation(x=x, fx=fx, w=w, I_c=I_c)
end
