using FiniteLineSource
using WGLMakie
using CairoMakie
using FastGaussQuadrature


points = [0.5, 1, 2, 3, 4]
# discretization = [20, 30, 30, 20]
discretization = [10, 10]

H = 10. /0.1
A = 0.
# line_limits=[0., 0.3, 0.5, 0.7, 1.]
line_limits=[0.,0.5, 1.]

function get_discretization(a,b,n)
    x, w = gausslegendre(n+1)   
    m = (b-a)/2
    c = (b+a)/2
    @. x = x * m + c 
end


get_discretization(line_limits[1]*H, line_limits[2]*H, discretization[1])

d = Vector{Float64}[]
for i=1:length(line_limits)-1
  x = get_discretization(line_limits[i]*H, line_limits[i+1]*H, discretization[i])
  push!(d, x)
end
d = vcat(d...)


##.
fig = Figure(size = (400,400))
# grid = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1])
lines!(ax, [0,0], [0,H]; color=:grey, linewidth = 2)
for l in line_limits*H
    lines!(ax, [-0.1,0.1], [l,l]; color=:indianred, linewidth = 2)
end
fig


scatter!(ax, zeros(length(d)), d, markersize = 15)
# xlims!(ax, -0,5, 0.5)
ax.aspect = .4
hidedecorations!(ax)

hidespines!(ax)
##.

text!(ax, 0.3, H; text = "0", fontsize = 20, align = (:center, :center))
text!(ax, 0.3, .5*H ; text = "$(.5*H)", fontsize = 20, align = (:center, :center))
text!(ax, 0.3, 0 ; text = "$(H)", fontsize = 20, align = (:center, :center))

text!(ax, -0.4, 0.5*H ; text = L"H/r_b", fontsize = 28, align = (:center, :center))

text!(ax, 0.4, 0.2*H ; text = L"N=11", fontsize = 22, align = (:center, :center))

xlims!(ax,-0.55,0.6)
fig


