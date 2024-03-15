using FiniteLineSource
using FiniteLineSource: convolve_step
using WGLMakie
using CairoMakie

function produce_plot(r, segment_points, segment_limits, label, color)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = PointToPoint(r=r)

    params = Constants(Δt = 3600., segment_limits=segment_limits, segment_points=segment_points)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, Δt = params.Δt, r = setup.r)

    abs_error = abs.(C-I) 
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)

    lines!(t, error, color=color)#, label=label, ylimits=(-20, 0))
end


rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (600,600))

segment_limits_1 = [0., 10.]
segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = "t (h)", ylabel = "log₁₀(ϵ)", title="r̃=$(Int(r/0.1))")
    produce_plot(r, [50], segment_limits_1, "n=50", :blue)    
    produce_plot(r, [100], segment_limits_1, "n=100", :red)       
    produce_plot(r, [200], segment_limits_1, "n=200", :green)       
    produce_plot(r, [400], segment_limits_1, "n=400", :orange)     
    produce_plot(r, [600], segment_limits_1, "n=600", :black)         
end
ff
