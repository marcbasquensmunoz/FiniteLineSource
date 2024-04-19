using FiniteLineSource
using FiniteLineSource: f_evolve_1!, f_evolve_2!
using WGLMakie
using CairoMakie
using Parameters
using LaTeXStrings

WGLMakie.activate!()

# compute big F the smooth function 

# compute_f stores the values for all time steps
function compute_f!(setup, Fx, q, x, fx, CΔt, params)
    n = length(q)
    for i = 1:n-1
        @inbounds f_evolve_1!(setup, fx, x, CΔt, q[i], params)
        @inbounds f_evolve_2!(setup, fx, x, q[i], params)
        @inbounds Fx[i+1,:] .= copy(fx)
    end
end

Δt = 3600.; α = 1e-6; kg = 3.; rb = 0.1; 
Δt̃ = Δt*α/rb^2
r=1.

setup = PointToPoint(r=r)
params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)
precomp = precompute_parameters(setup, params=params)
@unpack x, fx, v = precomp

guide(z) = exp(-10*z^2*Δt) * (1 - exp(-z^2*Δt)) / z
segments = adaptive_gk_segments(guide, 0., 10.)
segment_limits = sort(union([s.a for s in segments], [s.b for s in segments]))
#x, w = adaptive_gk(guide, 0., 10.)


CΔt = @. exp(-x^2*Δt̃)
q = [20*sin(2π*i/8760)+ 5 *sin(2π*i/8760)+ 5. for i=1:10000]
Fx = zeros(length(q),length(x))
@time compute_f!(setup, Fx, q, x, fx, CΔt, params)


##.
# Evolution of Fx in time 

ff = Figure(size = (600, 300))
ax = Axis(ff[1, 1], xlabel = latexstring(L"\zeta"), ylabel = latexstring(L"\bar{F}"), title=latexstring(L"\tilde{t} = 1"))

perm = sortperm(x)
lines!(ax, x[perm], Fx[10000,perm], color= :black, label = "1", linewidth=0.5)
scatter!(ax,  x[1:2:end], Fx[10000,:][1:2:end], color=:blue, markersize=4)
colors = [(i%2 == 0 ? :red : :grey, 0.07) for (i, lim) in enumerate(segment_limits[1:end-1])]
vspan!(segment_limits[1:end-1], segment_limits[2:end], color = colors)
xlims!(ax, 0, 10)
ax.xticks = 0:2:10
vlines!(segment_limits, colors= :blue, linewidth=1)

#vlines!(ax, segment_limits, color = :red)

ax2 = Axis(ff, bbox = BBox(327, 570, 130, 230), xticklabelsize = 12, yticklabelsize = 12, title = "zoomed view")
xlims!(ax2, 0, 0.5)
ax2.xticks = 0:0.1:0.5
lines!(ax2, x[perm], Fx[10000,perm], color= :black, strokewidth = 0.5, label = "1", linewidth=0.5)
scatter!(ax2,  x[1:2:end], Fx[10000,:][1:2:end], color=:blue, markersize=4)
vspan!(segment_limits[1:end-1], segment_limits[2:end], color = colors)
vlines!(segment_limits, colors= :blue, linewidth=1)

#Box(ff, color = :bisque, strokecolor = :blue, strokewidth = 1)

#axislegend(latexstring(L"\tilde{t}"); position= :rt, backgroundcolor = (:grey90, 0.25));

ff

CairoMakie.activate!()
save("figures/Fx_disc.pdf", ff)