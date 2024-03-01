using FiniteLineSource
using FiniteLineSource: convolve_fls_step, convolve_sls_step
using Plots


function produce_plot_fls(line_points, setup, label)
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    #setup = SegmentToPoint(D=10., H=20., z=15., r=1.)
    #setup = SegmentToSegment(D1=10., H1=20., D2=10., H2=20., r=1.)

    segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
    segment_points = 4 .* [10, 25, 10, 10, 10, 10]
    params = Constants(Δt = 3600., line_points=line_points, segment_limits=segment_limits, segment_points=segment_points)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_fls_step(q, Δt=params.Δt, r=setup.r, z=setup.z, D=setup.D, H=setup.H)
    #C = convolve_sls_step(q, Δt=params.Δt, r=setup.r, D1=setup.D1, H1=setup.H1, D2=setup.D2, H2=setup.H2)

    abs_error = abs.(C-I) 
    rel_error = abs.((C-I) ./ C) 

    plot!(1:length(q), log10.(abs_error), label=label, ylimits=(-20, 0))
end
#maximum(abs_error)
#plot(1:length(q), log10.(abs_error))

#maximum(rel_error)
#plot(1:length(q), log10.(rel_error))

rs = [0.1, 1, 10, 100]
plots = []
for r in rs
    f = plot(title="r̃=$(Int(r/0.1))", legend=r==100 ? true : false)
    produce_plot_fls(30, SegmentToPoint(D=10., H=20., z=15., r=r), "N=30")   
    produce_plot_fls(50, SegmentToPoint(D=10., H=20., z=15., r=r), "N=50")   
    produce_plot_fls(70, SegmentToPoint(D=10., H=20., z=15., r=r), "N=70")   
    produce_plot_fls(90, SegmentToPoint(D=10., H=20., z=15., r=r), "N=90")   
    produce_plot_fls(150, SegmentToPoint(D=10., H=20., z=15., r=r), "N=150")   

    xlabel!("t (h)")  
    ylabel!("log₁₀(ϵ)")
    push!(plots, f)
end
plot(plots..., layout=(2,2))
