using FiniteLineSource
using FiniteLineSource: convolve_step
using Plots

function produce_plot(r, segment_points, label)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = PointToPoint(r=r)

    segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
    #segment_points = 6 .* [10, 25, 10, 10, 10, 10]

    #segment_limits = [0., 10.]

    params = Constants(Δt = 3600., segment_limits=segment_limits, segment_points=segment_points)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, Δt = params.Δt, r = setup.r)

    abs_error = abs.(C-I) 
    rel_error = abs.((C-I) ./ C) 

    maximum(abs_error)
    plot!(1:length(q), log10.(abs_error), label=label, ylimits=(-20, 0))
    #plot(1:length(q), log10.(abs_error))

    #maximum(rel_error)
    #plot(1:length(q), log10.(rel_error))
end

rs = [0.1, 1, 10, 100]
plots = []
for r in rs
    f = plot(title="r̃=$(Int(r/0.1))", legend=r==100 ? true : false)
    produce_plot(r, [50], "n=50")    
    produce_plot(r, [100], "n=100")       
    produce_plot(r, [200], "n=200")       
    produce_plot(r, [400], "n=400")     
    produce_plot(r, [600], "n=600")         
    xlabel!("t (h)")  
    ylabel!("log₁₀(ϵ)")
    push!(plots, f)
end
plot(plots..., layout=(2,2))


rs = [0.1, 1, 10, 100]
plots = []
for r in rs
    f = plot(title="r̃=$(Int(r/0.1))", legend=r==100 ? true : false)
    produce_plot(r, [5, 20, 5, 5, 5, 5], "∑n=45")    
    produce_plot(r, 1 .* [10, 25, 10, 10, 10, 10], "∑n=75")       
    produce_plot(r, 2 .* [10, 25, 10, 10, 10, 10], "∑n=150")       
    produce_plot(r, 3 .* [10, 25, 10, 10, 10, 10], "∑n=225")     
    produce_plot(r, 4 .* [10, 25, 10, 10, 10, 10], "∑n=300")         
    xlabel!("t (h)")  
    ylabel!("log₁₀(ϵ)")
    push!(plots, f)
end
plot(plots..., layout=(2,2))
