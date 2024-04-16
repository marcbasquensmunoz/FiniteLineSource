using FiniteLineSource
using WGLMakie
using CairoMakie

function produce_plot_sts(r, label)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = SegmentToSegment(D1=5., H1=20., D2=5., H2=15., r=r)

    segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
    segment_points = 6 .* [10, 25, 10, 10, 10, 10]

    params = Constants(Δt = 3600., segment_limits=segment_limits, segment_points=segment_points)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, setup, params=params)

    abs_error = abs.(C-I) 
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)

    lines!(t ./ 8760, error, label=label)
end

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (800,600))

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel =  floor(Int, (i-1)/2) == 1 ? "t (years)" : "", ylabel = i%2 == 1 ? "log₁₀(ϵ)" : "", title="r̃=$(Int(r/0.1))", limits=(nothing, (-20, 0)))
    produce_plot_sts(r, "")    
end
ff

CairoMakie.activate!()
save("figures/error_sts.pdf", ff)