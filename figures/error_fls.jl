using FiniteLineSource
using WGLMakie
using CairoMakie

function produce_plot_fls(r, line_points)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = SegmentToPoint(D=5., H=50., r=r, z=30.)

    params = Constants(Δt = 3600., line_points=line_points, line_limits=[0., 0.3, 0.5, 0.7, 1.])

    precomp = precompute_parameters(setup, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, setup, params=params)

    abs_error = abs.(C-I) 
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)

    lines!(t ./ 8760, error, label="∑N=$(Int(sum(line_points)))")
end

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (800,600))

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel =  floor(Int, (i-1)/2) == 1 ? "t (years)" : "", ylabel = i%2 == 1 ? "log₁₀(ϵ)" : "", title="r̃=$(Int(r/0.1))", limits=(nothing, (-20, 0)))
    produce_plot_fls(r, 0.5 .* [20, 30, 30, 20])    
    produce_plot_fls(r, [20, 30, 30, 20])       
    produce_plot_fls(r, 2 .* [20, 30, 30, 20])       
    produce_plot_fls(r, 3 .* [20, 30, 30, 20])     
    produce_plot_fls(r, 4 .* [20, 30, 30, 20])       
    ff[2, 3] = Legend(ff, ax, "")
end
ff

CairoMakie.activate!()
save("figures/error_fls.pdf", ff)