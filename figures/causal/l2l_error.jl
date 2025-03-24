using FiniteLineSource
using FiniteLineSource: LineSource
using Makie
using CairoMakie

CairoMakie.activate!()

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

D = 0.
H = 150.

Nt = 1000

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)
q = [1. for t in 1:Nt]

function compute_error_Linf_l2l(;σ, ϵ, q, constants, D, H)
    Nt = length(q)
    sources = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=σ, y=0., D=D, H=H)]

    setup = SegmentToSegment(D1=D, H1=H, D2=D, H2=H, σ=σ)

    containers = FiniteLineSource.AsymptoticContainers(10)
    block = prepare_containers(setup, sources, ϵ, Nt, constants, containers);
    Ib = zeros(length(sources), Nt)
    evolve!(Ib, q, block)

    C = convolve_step(q, SegmentToSegmentOld(setup); params=constants)

    err = @. @views abs(Ib[1, :] - C)
    mean(err)
end


r_range = 0.7:0.5:10.
ϵ_range = [1e-4, 1e-6, 1e-8, 1e-10]#, 1e-12]

res = zeros(length(r_range), length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res[:, i] = [compute_error_Linf_l2l(σ=σσ, ϵ=ϵ, q=q, constants=constants, D=D, H=H) for σσ in r_range]
end

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"\sigma", ylabel = L"\log_{10} \Vert \epsilon  \Vert_{\infty}", title = L"\text{Error in the line to line case}; \ D_s = D_t = 0, \ H_s=H_t=150 m")

for (i, ϵ) in enumerate(ϵ_range)
    lines!(ax, r_range, log10.(res[:, i]), label = L"\epsilon = 10^{%$(Int(log10(ϵ)))}")
end
Makie.xlims!(ax, 0, 10)
Makie.ylims!(ax, -16., 0.)
axislegend(""; position= :rt, backgroundcolor = (:grey90, 0.25));

fig

save("figures/causal/l2l_error.pdf", fig)
