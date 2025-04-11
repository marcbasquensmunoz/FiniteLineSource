using FiniteLineSource
using FiniteLineSource: PointSource
using Makie
using CairoMakie

include("error_utils.jl")

function compute_error_Linf_p2p(;r, ϵ, q_gen, constants)
    sources = [PointSource(0., 0., 0.), PointSource(r, 0., 0.)]

    Nt = 40*8760
    q = q_gen.(1:Nt)
    block = prepare_containers(PointToPoint(r=r), sources, ϵ, Nt, constants, Q=maximum(q));
    Ib = zeros(length(sources), Nt)
    evolve!(Ib, q, block)

    setup = PointToPoint(r = r)
    C = convolve_step(q, setup; params=constants)

    err = @. @views abs(Ib[1, :] - C)
    maximum(err)
end

res_step = zeros(length(r̃_range), length(ϵ_range))
res_synth = zeros(length(r̃_range), length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res_step[:, i] = [compute_error_Linf_p2p(r=rr*rb, ϵ=ϵ, q_gen=q_step, constants=constants) for rr in r̃_range]
    @. res_synth[:, i] = [compute_error_Linf_p2p(r=rr*rb, ϵ=ϵ, q_gen=q_synth, constants=constants) for rr in r̃_range]
end

fig = create_error_plot(ϵ_range, r̃_range, res_step, res_synth; xlabel=:r, title= L"\text{Error in the point to point case}")

save("figures/causal/error/p2p_error.pdf", fig)
