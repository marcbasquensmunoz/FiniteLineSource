using FiniteLineSource
using FiniteLineSource: LineSource, AsymptoticContainers
using Makie
using CairoMakie

include("error_utils.jl")

function compute_error_Linf_l2p(;σ, ϵ, q, constants, D, H, z, containers)
    Nt = length(q)
    sources = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=σ, y=0., D=D, H=H)]

    setup = SegmentToPoint(D=D, H=H, z=z, σ=σ)

    block = prepare_containers(setup, sources, ϵ, Nt, constants, containers);
    Ib = zeros(length(sources), Nt)
    evolve!(Ib, q, block)

    C = convolve_step(q, setup; params=constants)

    err = @. @views abs(Ib[1, :] - C)
    if maximum(err) > ϵ
        @show σ, ϵ, D, H, z
    end
    maximum(err)
end

D = 0.
H = 150.
z = D + H/2

res_step = zeros(length(r_range), length(ϵ_range))
res_synth = zeros(length(r_range), length(ϵ_range))
containers = AsymptoticContainers(10)

for (i, ϵ) in enumerate(ϵ_range)
    @info "Computing errors for ϵ=$ϵ"
    @. res_step[:, i] = [compute_error_Linf_l2p(σ=σσ, ϵ=ϵ, q=q_step, constants=constants, D=D, H=H, z=z, containers=containers) for σσ in r_range]
    #@. res_synth[:, i] = [compute_error_Linf_l2p(σ=σσ, ϵ=ϵ, q=q_synth, constants=constants, D=D, H=H, z=z, containers=containers) for σσ in r_range]
end

fig = create_error_plot(ϵ_range, r_range, res_step, res_synth; rb=rb, title= L"\text{Error in the line to point case}; \ \tilde{D} = %$(Int(D/rb)), \ \tilde{H}=%$(Int(H/rb))")

save("figures/causal/error/l2p_error.pdf", fig)
