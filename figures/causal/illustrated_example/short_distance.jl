using FiniteLineSource
using FiniteLineSource: LineSource, LineIntegralParams
using Makie
using CairoMakie
using Colors

######################
# Parameters 
######################
Nt = 8760

ϵ = 1e-6

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

D = 0.
H = 150.
σ = 5.
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

################################################
# Make plots of individual blocks (except first)
################################################

function get_block(block, B)
    ζ = 0:0.0001:(B == 1 ? 15. : maximum(block.ζ[block.ranges[B-1]]))
    block_integrand = zeros(length(ζ))
    block_integrand_full = zeros(length(ζ))
    NB = vcat(0, block.N)
    s = [sum([q[i] * exp(-ζζ^2 * constants.Δt̃ * (Nt - i)) for i in Nt-NB[B+1]+1:Nt-NB[B]]) for ζζ in ζ]
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
        lineparams_full = LineIntegralParams(setup, ζζ/rb, ϵ)
        block_integrand[k] = (1 - exp(-ζζ^2 * constants.Δt̃)) / ζζ * FiniteLineSource.compute_line_integral(lineparams; containers=containers) * s[k]
        block_integrand_full[k] = (1 - exp(-ζζ^2 * constants.Δt̃)) / ζζ * FiniteLineSource.compute_line_integral(lineparams_full; containers=containers) * s[k]
    end
    ζ, block_integrand, block_integrand_full, 2h
end

ζ = [zeros(0) for _ in 1:length(block.N)]
I_region = [zeros(0) for _ in 1:length(block.N)]
I_full = [zeros(0) for _ in 1:length(block.N)]
H_eff = zeros(length(block.N))

for i in 1:length(block.N)
    ζ[i], I_region[i], I_full[i], H_eff[i] = get_block(block, i)
end

primary_colors = Makie.wong_colors()
lighten(color) = RGBA(color.r, color.g, color.b, 0.5)
secondary_colors = lighten.(primary_colors)

block_figs = Figure[]

for i in 2:length(block.N)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\zeta", title = "Integrand for block $i")

    lines!(ax, ζ[i], I_region[i], label = "Reduced integration region", color=primary_colors[i-1])
    lines!(ax, ζ[i], I_full[i], label = "Full source integration", color=secondary_colors[i-1])

    axislegend("", position = :rt)

    push!(block_figs, fig)
    save("$(@__DIR__)/block_$i.pdf", fig)
end


############################################
# Make plot of first block 
############################################


fig_fist = Figure()
ax = Axis(fig_fist[1, 1], xlabel = L"\zeta", title = "Integrand for block 1")
lines!(ax, ζ[1], I_region[1], label = "Reduced integration region")

save("$(@__DIR__)/first_block.pdf", fig_fist)



###############################################
# Make image of lines with integration regions 
###############################################

fig_lines = Figure()
x_pos = [75*(i-1) for i in 1:length(block.N)-1]
block_height_margin = 3.
line_limit_width = 1.
text_margin = 10.
figure_margin = 20.

ax = Axis(fig_lines[1, 1])
hidedecorations!(ax)
hidespines!(ax) 

for (i, x) in enumerate(x_pos)
    B = i+1
    # Draw dource line
    lines!(ax, [x, x], [-H/2,H/2], color=:black)
    lines!(ax, [x-line_limit_width, x+line_limit_width], [-H/2,-H/2], color=:black)
    lines!(ax, [x-line_limit_width, x+line_limit_width], [H/2,H/2], color=:black)
    
    # Draw target point
    scatter!(ax, [x+σ], [0], color=:black)

    # Draw reduced integration region
    HB = H_eff[B]
    lines!(ax, [x, x], [-HB/2,HB/2], color=primary_colors[i], linewidth=6)
    text!(ax, x-block_height_margin, 0, text=L"%$(round(HB, digits=1))m", align = (:right, :center), fontsize=20, color=primary_colors[i])

    text!(ax, x, H/2+text_margin, text="Block $B", align = (:center, :center), fonts =(; regular= "arial"), fontsize=20)    
end
text!(ax, sum(x_pos)/2, -H/2-text_margin, text=L"H=%$(Int(round(H, digits=0)))m, \ σ=%$(Int(round(σ, digits=0)))m ", align = (:center, :center), fonts =(; regular= "arial"), fontsize=20)
    
xlims!(ax, x_pos[1]-figure_margin, x_pos[end]+figure_margin)
fig_lines   
save("$(@__DIR__)/line_sources_blocks.pdf", fig_lines)