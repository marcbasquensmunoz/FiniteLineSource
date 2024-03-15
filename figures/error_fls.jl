using FiniteLineSource
using FiniteLineSource: convolve_fls_step
using WGLMakie
using CairoMakie

function produce_plot_fls(r, segment_points, label)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = SegmentToPoint(D=5., H=50., r=r, z=30.)

    segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]

    params = Constants(Δt = 3600., segment_limits=segment_limits, segment_points=segment_points)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_fls_step(q, Δt=3600., r=setup.r, z=setup.z, D=setup.D, H=setup.H, α=params.α, kg=params.kg)

    abs_error = abs.(C-I) 
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)

    lines!(t, error, label=label)
end

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (600,600))

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = "t (h)", ylabel = "log₁₀(ϵ)", title="r̃=$(Int(r/0.1))")
    produce_plot_fls(r, [5, 20, 5, 5, 5, 5], "∑n=45")    
    produce_plot_fls(r, 1 .* [10, 25, 10, 10, 10, 10], "∑n=75")       
    produce_plot_fls(r, 2 .* [10, 25, 10, 10, 10, 10], "∑n=150")       
    produce_plot_fls(r, 3 .* [10, 25, 10, 10, 10, 10], "∑n=225")     
    produce_plot_fls(r, 4 .* [10, 25, 10, 10, 10, 10], "∑n=300")       
end
ff