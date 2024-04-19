using FiniteLineSource
using FiniteLineSource: convolve_step
using WGLMakie
using CairoMakie

function produce_plot(r, n)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = PointToPoint(r=r)

    params = Constants(Δt = 3600.)

    n_tot = [0]
    precomp = precompute_parameters(setup, params=params, n=n, n_tot=n_tot)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, Δt = params.Δt, r = setup.r)

    abs_error = abs.(C-I) 
    t = 1:length(q)
    abs_error[abs_error .== 0] .= 10^-32
    error = log10.(abs_error)

    lines!(t ./ 8760, error, label="∑n=$(n_tot[1])")
end

####

rs = [0.1, 1, 10, 100]
plots = []
ff = Figure(size = (800,500))

for (i, r) in enumerate(rs)
    ax = Axis(ff[floor(Int, (i-1)/2) + 1, (i-1)%2+1], xlabel = floor(Int, (i-1)/2) == 1 ? "t (years)" : "", ylabel = i%2 == 1 ? "log₁₀(ϵ)" : "", title="r̃=$(Int(r/0.1))", limits=(nothing, (-20, 0)))

    produce_plot(r, 4)    
    produce_plot(r, 8)       
    produce_plot(r, 12)       
    produce_plot(r, 16)     
    produce_plot(r, 20)      
    ff[2, 3] = Legend(ff, ax, "")   
end
ff

CairoMakie.activate!()
save("figures/error_ps.pdf", ff)
