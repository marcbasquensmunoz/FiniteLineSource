using FiniteLineSource
using FiniteLineSource: LineSource
using Makie
using CairoMakie

include("error_utils.jl")

function compute_error_Linf_l2l(;σ, ϵ, q_gen, constants, D, H, containers)
    sources = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=σ, y=0., D=D, H=H)]
    setup = SegmentToSegment(D1=D, H1=H, D2=D, H2=H, σ=σ)

    Nt = 8760*15
    q = q_gen.(1:Nt)
    block = prepare_containers(setup, sources, ϵ, Nt, constants, containers, Q=maximum(q));
    Ib = zeros(length(sources), Nt)
    evolve!(Ib, q, block)

    C = convolve_step(q, SegmentToSegmentOld(setup); params=constants, ϵ=ϵ/10)

    err = abs.( (Ib[1, :] - C))
    if maximum(err) > ϵ
        @show σ, ϵ, D, H
    end
    maximum(err)
end

D = 0.
H = 150.

res_step = zeros(length(r̃_range), length(ϵ_range))
res_synth = zeros(length(r̃_range), length(ϵ_range))

containers = AsymptoticContainers(10)

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res_step[:, i] = [compute_error_Linf_l2l(σ=σσ*rb, ϵ=ϵ, q_gen=q_step, constants=constants, D=D, H=H, containers=containers) for σσ in r̃_range]
    @. res_synth[:, i] = [compute_error_Linf_l2l(σ=σσ*rb, ϵ=ϵ, q_gen=q_synth, constants=constants, D=D, H=H, containers=containers) for σσ in r̃_range]
end

fig = create_error_plot(ϵ_range, r̃_range, res_step, res_synth; title= L"\text{Error in the line to line case}; \ \tilde{D}_s = \tilde{D}_t = %$(Int(D/rb)), \ \tilde{H}_s=\tilde{H}_t=%$(Int(H/rb))")

save("figures/causal/error/l2l_error.pdf", fig)
