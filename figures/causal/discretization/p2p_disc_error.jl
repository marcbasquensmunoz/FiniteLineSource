using FiniteLineSource
using FiniteLineSource: LineSource
using Makie
using CairoMakie
using Parameters

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

r = 1.
setup = PointToPoint(r=r)

Nt = 10000

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

function count_ζ_original(;setup, ϵ, constants, Nt)
    precomp = precompute_parameters(setup, params=constants, ϵ=ϵ);
    return length(precomp.x), ϵ
end
function count_ζ_block(;setup, ϵ, constants, Nt)
    @unpack r = setup
    sources = [(x=0., y=0., z=0.), (x=r, y=0., z=0.),]
    block = prepare_containers(setup, sources, ϵ, Nt, constants);
    return length(block.ζ), ϵ
end

ϵ_range = 10. .^ collect(0:-1:-13)

points_original = zeros(2, length(ϵ_range))
points_block = zeros(2, length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Counting ζ points for ϵ=$ϵ"
    points_block[:, i] .= count_ζ_block(setup=setup, ϵ=ϵ, constants=constants, Nt=Nt)
    points_original[:, i] .= count_ζ_original(setup=setup, ϵ=ϵ, constants=constants, Nt=Nt)
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"\log_{10} \ N_D", ylabel = L"\log_{10} \ \Vert \epsilon  \Vert_{\infty}", title = L"\text{Discretization: point to point}")

lines!(ax, log10.(points_block[1, :]), log10.(points_block[2, :]), label = "Blocks method")
lines!(ax, log10.(points_original[1, :]), log10.(points_original[2, :]), label = "Original method")

Makie.xlims!(ax, 1, 4)
Makie.ylims!(ax, -16, 0)
axislegend(""; position= :rt, backgroundcolor = (:grey90, 0.25));

fig

save("$(@__DIR__)/p2p_disc.pdf", fig)
