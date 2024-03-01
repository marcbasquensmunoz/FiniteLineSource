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
q = [20*sin(2π*i/8760)+ 5 *sin(2π*i/8760)+ 5. for i=1:200000]
Fx = zeros(length(q),length(x))
@time compute_f!(Fx,q,x,fx,CΔt)

##.
# Evolution of Fx in time 

ff = Figure(size = (600,250))
ax = Axis(ff[1, 1], xlabel = latexstring(L"\zeta"), ylabel = latexstring(L"\bar{F}"))

lines!(ax, x, Fx[100,:], color= :green, label = "0.01")
lines!(ax, x, Fx[1000,:], color= :red, label = "0.1")
lines!(ax, x, Fx[10000,:], color= :blue, label = "1")
lines!(ax, x, Fx[100000,:], color= :dodgerblue, label = "10")
xlims!(ax, 0, .5)
axislegend(latexstring(L"\tilde{t}"); position= :rt, backgroundcolor = (:grey90, 0.25));

ff

CairoMakie.activate!()
save("figures/Fx.pdf", ff)