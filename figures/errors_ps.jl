using FiniteLineSource
using FiniteLineSource: convolve_step
using WGLMakie
using CairoMakie

function produce_plot(r, segment_points, segment_limits, label)        
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

    lines!(t ./ 8760, error, label=label)
end

####
## Plots with m=1

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (800,500))

segment_limits = [0., 10.]

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = floor(Int, (i-1)/2) == 1 ? "t (years)" : "", ylabel = i%2 == 1 ? "log₁₀(ϵ)" : "", title="r̃=$(Int(r/0.1))", limits=(nothing, (-20, 0)))

    produce_plot(r, [50], segment_limits, "n=50")    
    produce_plot(r, [100], segment_limits, "n=100")       
    produce_plot(r, [200], segment_limits, "n=200")       
    produce_plot(r, [400], segment_limits, "n=400")     
    produce_plot(r, [600], segment_limits, "n=600")      
    ff[2, 3] = Legend(ff, ax, "")   
end
ff

CairoMakie.activate!()
save("figures/error_m1.pdf", ff)
##



####
## Plots with m=6

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (800,500))

segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = floor(Int, (i-1)/2) == 1 ? "t (years)" : "", ylabel = i%2 == 1 ? "log₁₀(ϵ)" : "", title="r̃=$(Int(r/0.1))", limits=(nothing, (-20, 0)))
    produce_plot(r, [10, 15, 5, 5, 5, 5], segment_limits, "∑n=45")    
    produce_plot(r, [10, 25, 10, 10, 10, 10], segment_limits, "∑n = 75")       
    produce_plot(r, 2 .* [10, 25, 10, 10, 10, 10], segment_limits, "∑n = 150")       
    produce_plot(r, 3 .* [10, 25, 10, 10, 10, 10], segment_limits, "∑n = 225")     
    produce_plot(r, 4 .* [10, 25, 10, 10, 10, 10], segment_limits, "∑n = 300")   
    ff[2, 3] = Legend(ff, ax, "")    
end
ff

CairoMakie.activate!()
save("figures/error_m6.pdf", ff)
##