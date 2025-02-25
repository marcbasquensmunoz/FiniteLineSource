using FiniteLineSource
using Makie
using CairoMakie

CairoMakie.activate!()

ϵ = 1e-12

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

Nt = 100

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)
q = [1. for t in 1:Nt]

function compute_error_Linf(;r, ϵ, q, constants)
    Nt = length(q)
    sources = [(0., 0., 0.), (r, 0., 0.)]

    block = prepare_containers(PointToPoint(r=r), sources, ϵ, Nt, constants);
    Ib = zeros(length(sources), Nt)
    evolve!(Ib, q, block)

    setup = PointToPoint(r = r)
    C = convolve_step(q, setup; params=constants)

    err = @. @views abs(Ib[1, :] - C)
    maximum(err)
end


r_range = 0.7:0.01:10.
ϵ_range = [1e-6, 1e-8, 1e-10, 1e-12, 1e-14]

res = zeros(length(r_range), length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res[:, i] = [compute_error_Linf(r=rr, ϵ=ϵ, q=q, constants=constants) for rr in r_range]
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"r", ylabel = L"\log_{10} \Vert \epsilon _{\infty} \Vert", title = L"\text{Error in the point to point case}")

for (i, ϵ) in enumerate(ϵ_range)
    lines!(ax, r_range, log10.(res[:, i]), label = L"\epsilon = 10^{%$(Int(log10(ϵ)))}")
end
xlims!(ax, 0, 10)
ylims!(ax, -16., 0.)
axislegend(""; position= :rt, backgroundcolor = (:grey90, 0.25));

fig

save("figures/causal/p2p_error.pdf", fig)
