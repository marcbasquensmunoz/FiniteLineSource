using FiniteLineSource
using FiniteLineSource: f_evolve_1!, f_evolve_2!
using WGLMakie
using CairoMakie
using Parameters
using LaTeXStrings

WGLMakie.activate!()

# compute big F the smooth function 

# compute_f stores the values for all time steps
function compute_f!(Fx, q, x, fx, CΔt)
    n = length(q)
    for i = 1:n-1
        @inbounds f_evolve_1!(fx, x, CΔt, q[i])
        @inbounds f_evolve_2!(fx, x, q[i])
        @inbounds Fx[i+1,:] .= copy(fx)
    end
end

Δt = 3600.; α = 1e-6; kg = 3.; rb = 0.1; 
Δt̃ = Δt*α/rb^2
r=1.

setup = PointToPoint(r=r)
segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
segment_points = 6 .* [10, 25, 10, 10, 10, 10]
params = Constants(Δt=Δt, α=α, kg=kg, rb=rb, segment_limits=segment_limits, segment_points=segment_points)
prealloc = Preallocation(setup, params) 

precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
@unpack x,fx,v = precomp

CΔt = @. exp(-x^2*Δt̃)
q = [20*sin(2π*i/8760)+ 5 *sin(2π*i/8760)+ 5. for i=1:10000]
Fx = zeros(length(q),length(x))
@time compute_f!(Fx,q,x,fx,CΔt)

##.
# Evolution of Fx in time 

ff = Figure(size = (600, 300))
ax = Axis(ff[1, 1], xlabel = latexstring(L"\zeta"), ylabel = latexstring(L"\bar{F}"), title=latexstring(L"\tilde{t} = 1"))

lines!(ax, x, Fx[10000,:], color= :black, label = "1", linewidth=0.5)
scatter!(ax,  x[1:2:end], Fx[10000,:][1:2:end], color=:blue, markersize=4)
colors = [(i%2 == 0 ? :grey : :red, 0.07) for (i, lim) in enumerate(segment_limits[1:end-1])]
vspan!(segment_limits[1:end-1], segment_limits[2:end], color = colors)
xlims!(ax, 0, 10)
ax.xticks = 0:2:10
vlines!(segment_limits, colors= :blue, linewidth=1)

#vlines!(ax, segment_limits, color = :red)

ax2 = Axis(ff, bbox = BBox(320, 570, 130, 230), xticklabelsize = 12, yticklabelsize = 12, title = "zoomed view")
xlims!(ax2, 0, 0.5)
ax2.xticks = 0:0.1:0.5
lines!(ax2, x, Fx[10000,:], color= :black, strokewidth = 0.5, label = "1", linewidth=0.5)
scatter!(ax2,  x[1:2:end], Fx[10000,:][1:2:end], color=:blue, markersize=4)
vspan!(segment_limits[1:end-1], segment_limits[2:end], color = colors)
vlines!(segment_limits, colors= :blue, linewidth=1)

#Box(ff, color = :bisque, strokecolor = :blue, strokewidth = 1)

#axislegend(latexstring(L"\tilde{t}"); position= :rt, backgroundcolor = (:grey90, 0.25));

ff

CairoMakie.activate!()
save("figures/Fx_disc.pdf", ff)