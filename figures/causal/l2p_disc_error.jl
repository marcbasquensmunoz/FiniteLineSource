using FiniteLineSource
using FiniteLineSource: LineSource
using Makie
using CairoMakie
using Parameters

CairoMakie.activate!()

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

D = 0.
H = 150.
z = 75.
σ = 1.

Nt = 10000
setup = SegmentToPoint(D=D, H=H, z=z, σ=σ)

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

function count_ζ_original(;setup, n, constants, Nt)
    params = Constants(Δt=constants.Δt, α=constants.α, kg=constants.kg, rb=constants.rb, line_limits=[0., 0.3, 0.6, 1.], line_points=n .* [1, 1, 1])
    precomp = precompute_parameters(setup, params=params);
    I = zeros(Nt)
    q = ones(Nt)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
    C = convolve_step(q, setup; params=params)
    err = abs.(C .- I)
    return length(precomp.x) * sum(params.line_points), maximum(err)
end
function count_ζ_block(;setup, ϵ, constants, Nt)
    @unpack D, H, σ = setup
    sources = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=σ, y=0., D=D, H=H)]
    Nt = FiniteLineSource.compute_N(sqrt(σ^2 + (D-z)^2), ϵ, constants, Nt)
    containers = AsymptoticContainers(10)
    block = prepare_containers(setup, sources, ϵ, Nt, constants, containers);
    return length(block.ζ), ϵ
end

ϵ_range = 10. .^ collect(0:-1:-14)
n_range = Int.(floor.(10. .^ collect(0.:1/3:3)))

points_original = zeros(2, length(n_range))
points_block = zeros(2, length(ϵ_range))

for (i, ϵ) in enumerate(ϵ_range)
    @info "Counting ζ points in the block method for ϵ=$ϵ"
    points_block[:, i] .= count_ζ_block(setup=setup, ϵ=ϵ, constants=constants, Nt=Nt)
end
for (j, nn) in enumerate(n_range)
    @info "Counting ζ points in the original method for n=$nn"
    points_original[:, j] .= count_ζ_original(setup=setup, n=nn, constants=constants, Nt=Nt)
end


fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"\log_{10} \ N_D", ylabel = L"\log_{10} \ \Vert \epsilon  \Vert_{\infty}", title = L"\text{Discretization: line to point}")

lines!(ax, log10.(points_block[1, :]), log10.(points_block[2, :]), label = "Blocks method")
lines!(ax, log10.(points_original[1, :]), log10.(points_original[2, :]), label = "Original method")

Makie.xlims!(ax, 1, 7)
Makie.ylims!(ax, -16, 0)
axislegend(""; position= :rt, backgroundcolor = (:grey90, 0.25));

fig

save("figures/causal/l2p_disc.pdf", fig)
