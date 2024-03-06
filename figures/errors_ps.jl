using FiniteLineSource
using FiniteLineSource: convolve_step
using WGLMakie
using CairoMakie
using Polynomials

function produce_plot(r, segment_points, label, color)        
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
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)
    fit_f = fit(error, t, 21)

    lines!(t, fit_f.(t), color=color)#, label=label, ylimits=(-20, 0))
end


rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (600,600))

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = "t (h)", ylabel = "log₁₀(ϵ)", title="r̃=$(Int(r/0.1))")
    produce_plot(r, [50], "n=50", :blue)    
    #produce_plot(r, [100], "n=100", :red)       
    #produce_plot(r, [200], "n=200", :green)       
    #produce_plot(r, [400], "n=400", :orange)     
    #produce_plot(r, [600], "n=600", :black)         
end
ff


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
