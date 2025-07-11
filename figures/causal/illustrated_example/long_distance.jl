using FiniteLineSource
using FiniteLineSource: LineSource, LineIntegralParams
using Makie
using CairoMakie
using Colors
using Roots
using SpecialFunctions

######################
# Parameters 
######################
Nt = 8760*20

ϵ = 1e-6

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

D = 0.
H = 150.
σ = 50.
setup = SegmentToPoint(D=D, H=H, z=D+H/2, σ=σ)

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = ones(Nt)

sources = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=σ, y=0., D=D, H=H)]


############################################
# Simulate with block methods 
############################################
containers = AsymptoticContainers(10)
block = prepare_containers(setup, sources, ϵ, Nt, constants, containers, Q=maximum(q));
Ib = zeros(length(sources), Nt)
evolve!(Ib, q, block)

############################################
# Simulate with block methods 
############################################
function get_block(block, B)
    ζmax_block = (B == 1 ? 15. : maximum(block.ζ[block.ranges[B-1]]))
    ζ = 0:0.0001:ζmax_block
    ζmax = find_zero(b -> rb/σ / b * sqrt(π/constants.Δt̃) * erfc(b*sqrt(constants.Δt̃)) - ϵ, 1.)
    ζ_full = 0:0.0001:ζmax
    block_integrand = zeros(length(ζ))
    block_integrand_full = zeros(length(ζ_full))
    NB = vcat(0, block.N)
    #s = [sum([q[i] * exp(-ζζ^2 * constants.Δt̃ * (Nt - i)) for i in Nt-NB[B+1]+1:Nt-NB[B]]) for ζζ in ζ]
    #full_s = [sum([q[i] * exp(-ζζ^2 * constants.Δt̃ * (Nt - i)) for i in 1:Nt]) for ζζ in ζ_full]
    # Cheating since q = 1 for faster evaluation
    s = [ (exp(-ζζ^2 * constants.Δt̃ * NB[B]) - exp(-ζζ^2 * constants.Δt̃ * (NB[B+1]+1))) / (1 - exp(-ζζ^2 * constants.Δt̃)) for ζζ in ζ]
    full_s = [ (exp(-ζζ^2 * constants.Δt̃) - exp(-ζζ^2 * constants.Δt̃ * Nt)) / (1 - exp(-ζζ^2 * constants.Δt̃)) for ζζ in ζ_full]
    z_eval = D + H/2
    edge = max(abs(z_eval-D), abs(z_eval-D-H))
    local maxh
    try
        maxh = find_zero(h -> FiniteLineSource.I_L2P(h, NB[B+1], setup, constants, ϵ) - ϵ, σ)
    catch
        maxh = edge
    end
    ND = sqrt(σ^2 + min(edge^2, maxh^2))
    h = sqrt(ND^2 - σ^2)
    D_eval = max(D, z_eval - h)
    H_eval = min(D + H, z_eval + h) - D_eval
    setup_cut = SegmentToPoint(D=D_eval, H=H_eval, z=z_eval, σ=σ)
    for (k, ζζ) in enumerate(ζ)
        lineparams = LineIntegralParams(setup_cut, ζζ/rb, ϵ)
        block_integrand[k] = (1 - exp(-ζζ^2 * constants.Δt̃)) / ζζ * FiniteLineSource.compute_line_integral(lineparams; containers=containers) * s[k]
    end
    for (k, ζζ) in enumerate(ζ_full)
        lineparams = LineIntegralParams(setup_cut, ζζ/rb, ϵ)
        block_integrand_full[k] = (1 - exp(-ζζ^2 * constants.Δt̃)) / ζζ * FiniteLineSource.compute_line_integral(lineparams; containers=containers) * full_s[k]
    end
    ζ, block_integrand, ζ_full, block_integrand_full
end

ζ = [zeros(0) for _ in 1:length(block.N)]
ζ_full = [zeros(0) for _ in 1:length(block.N)]
I_block = [zeros(0) for _ in 1:length(block.N)]
I_full = [zeros(0) for _ in 1:length(block.N)]

for i in 2:length(block.N)
    ζ[i], I_block[i], ζ_full[i], I_full[i] = get_block(block, i)
end

primary_colors = Makie.wong_colors()
lighten(color) = RGBA(color.r, color.g, color.b, 0.5)
secondary_colors = lighten.(primary_colors)

block_figs = Figure[]

for i in 2:length(block.N)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\zeta", title = "Integrand")

    ζ_range = 1:div(length(ζ_full[i]),4)

    lines!(ax, ζ[i], I_block[i], label = "Loads in block $i", color=primary_colors[i-1])
    lines!(ax, ζ_full[i][ζ_range], I_full[i][ζ_range], label = "Full load history", color=secondary_colors[i-1])

    axislegend("", position = :rb)

    ax_inset = Axis(fig[1, 1],
        width=Relative(0.5),
        height=Relative(0.3),
        halign=0.9,
        valign=0.9,
        title="Zoomed View")
    
    n = findfirst(x -> x > maximum(ζ[i]), ζ_full[i])
    lines!(ax_inset, ζ[i], I_block[i], color=primary_colors[i-1])
    lines!(ax_inset, ζ_full[i][1:n], I_full[i][1:n], color=secondary_colors[i-1])

    push!(block_figs, fig)
    save("$(@__DIR__)/long_block_$i.pdf", fig)
end


lines( ζ_full[2][1:50], I_block[2][1:50])