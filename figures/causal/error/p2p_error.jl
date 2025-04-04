using FiniteLineSource
using Makie
using CairoMakie

include("error_utils.jl")

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

res_step = zeros(length(r_range), length(ϵ_range))
res_synth = zeros(length(r_range), length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res_step[:, i] = [compute_error_Linf(r=rr, ϵ=ϵ, q=q_step, constants=constants) for rr in r_range]
    #@. res_synth[:, i] = [compute_error_Linf(r=rr, ϵ=ϵ, q=q_synth, constants=constants) for rr in r_range]
end

fig = create_error_plot(ϵ_range, r_range, res_step, res_synth; xlabel=:r, rb=rb, title= L"\text{Error in the point to point case}")

save("figures/causal/error/p2p_error.pdf", fig)
