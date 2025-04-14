using FiniteLineSource
using FiniteLineSource: PointSource, LineSource
using BenchmarkTools
using Interpolations

Nt = 8760*20

Δt = 3600.
α = 1e-6
kg = 3.
rb = 0.1

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [1. for t in 1:Nt]
ϵ_range = 10. .^ collect(-2:-2:-12)
n_range = Int.(floor.(10. .^ collect(0.:1/3:2+2/3)))

B = 1.
H = 150.
D = 0.
point_positions = [PointSource(0., 0., 0.), PointSource(B, 0., 0.)]
line_positions = [LineSource(x=0., y=0., D=D, H=H), LineSource(x=B, y=0., D=D, H=H)]

setup_p2p = PointToPoint(r = B)
setup_l2p = SegmentToPoint(D=D, H=H, z=D+H/2, σ=B)
setup_l2l = SegmentToSegment(D1=D, H1=H, D2=D, H2=H, σ=B)

precomputation = zeros(length(ϵ_range), 6)
simulation = zeros(length(ϵ_range), 6)
containers = AsymptoticContainers(10)

Ib = zeros(2, Nt)
Inh = zeros(Nt)

function interpolate_error(n_range, constants, Nt, setup) 
    n_eff = zeros(length(n_range))
    err = zeros(length(n_range))

    for (i, n) in enumerate(n_range)
        @info "Computing errors for n=$n"
        params = Constants(Δt=constants.Δt, α=constants.α, kg=constants.kg, rb=constants.rb, line_limits=[0., 0.3, 0.6, 1.], line_points=n .* [1, 1, 1])
        precomp = precompute_parameters(setup, params=params);
        I = zeros(Nt)
        q = ones(Nt)
        compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
        C = convolve_step(q, setup; params=params)
        err[i] = maximum(abs.(C .- I))
        n_eff[i] = sum(params.line_points)
    end
    perm = sortperm(err)
    LinearInterpolation(err[perm], n_eff[perm])
end

N_l2p = interpolate_error(n_range, constants, Nt, setup_l2p)
N_l2l = interpolate_error(n_range, constants, Nt, SegmentToSegmentOld(setup_l2l))

for (i, ϵ) in enumerate(ϵ_range)
    @show ϵ

    # Point to point
    b_block = @benchmark prepare_containers($setup_p2p, $point_positions, $ϵ, $Nt, $constants);
    block = prepare_containers(setup_p2p, point_positions, ϵ, Nt, constants);
    b_evolve = @benchmark evolve!($Ib, $q, $block)

    b_precomp = @benchmark precompute_parameters($setup_p2p, params=$constants, ϵ=$ϵ);
    precomp = precompute_parameters(setup_p2p, params=constants, ϵ=ϵ);
    b_compute = @benchmark compute_integral_throught_history!($setup_p2p, I=$Inh, q=$q, precomp=$precomp, params=$constants)

    precomputation[i, 1] = minimum(b_precomp.times)
    simulation[i, 1] = minimum(b_compute.times)
    precomputation[i, 2] = minimum(b_block.times)
    simulation[i, 2] = minimum(b_evolve.times)

    # Line to point
    b_block = @benchmark prepare_containers($setup_l2p, $line_positions, $ϵ, $Nt, $constants, $containers);
    block = prepare_containers(setup_l2p, line_positions, ϵ, Nt, constants, containers);
    b_evolve = @benchmark evolve!($Ib, $q, $block)

    l2p_constants = Constants(Δt=constants.Δt, α=constants.α, kg=constants.kg, rb=constants.rb, line_limits=[0., 0.3, 0.6, 1.], line_points=Int(ceil(N_l2p(ϵ))) .* [1, 1, 1])
    b_precomp = @benchmark precompute_parameters($setup_l2p, params=$l2p_constants, ϵ=$ϵ);
    precomp = precompute_parameters(setup_l2p, params=l2p_constants, ϵ=ϵ);
    b_compute = @benchmark compute_integral_throught_history!($setup_l2p, I=$Inh, q=$q, precomp=$precomp, params=$l2p_constants)

    precomputation[i, 3] = minimum(b_precomp.times)
    simulation[i, 3] = minimum(b_compute.times)
    precomputation[i, 4] = minimum(b_block.times)
    simulation[i, 4] = minimum(b_evolve.times)

    # Line to line
    b_block = @benchmark prepare_containers($setup_l2l, $line_positions, $ϵ, $Nt, $constants, $containers);
    block = prepare_containers(setup_l2l, line_positions, ϵ, Nt, constants, containers);
    b_evolve = @benchmark evolve!($Ib, $q, $block)

    l2l_constants = Constants(Δt=constants.Δt, α=constants.α, kg=constants.kg, rb=constants.rb, line_limits=[0., 0.3, 0.6, 1.], line_points=Int(ceil(N_l2l(ϵ))) .* [1, 1, 1])
    b_precomp = @benchmark precompute_parameters(SegmentToSegmentOld($setup_l2l), params=$l2l_constants, ϵ=$ϵ);
    precomp = precompute_parameters(SegmentToSegmentOld(setup_l2l), params=l2l_constants, ϵ=ϵ);
    b_compute = @benchmark compute_integral_throught_history!(SegmentToSegmentOld($setup_l2l), I=$Inh, q=$q, precomp=$precomp, params=$l2l_constants)

    precomputation[i, 5] = minimum(b_precomp.times)
    simulation[i, 5] = minimum(b_compute.times)
    precomputation[i, 6] = minimum(b_block.times)
    simulation[i, 6] = minimum(b_evolve.times)
end


simulation_no_update = zeros(length(ϵ_range), 3)

for (i, ϵ) in enumerate(ϵ_range)
    # For this computation, comment out the lines corresponding to the update step in the evolve! function
    @show ϵ

    # Point to point
    block = prepare_containers(setup_p2p, point_positions, ϵ, Nt, constants);
    b_evolve = @benchmark evolve!($Ib, $q, $block)
    simulation_no_update[i, 1] = minimum(b_evolve.times)

    # Line to point
    block = prepare_containers(setup_l2p, line_positions, ϵ, Nt, constants, containers);
    b_evolve = @benchmark evolve!($Ib, $q, $block)
    simulation_no_update[i, 2] = minimum(b_evolve.times)

    # Line to line
    block = prepare_containers(setup_l2l, line_positions, ϵ, Nt, constants, containers);
    b_evolve = @benchmark evolve!($Ib, $q, $block)
    simulation_no_update[i, 3] = minimum(b_evolve.times)
end

convolution = zeros(length(ϵ_range), 3)

for (i, ϵ) in enumerate(ϵ_range)
    @show ϵ

    # Point to point
    b_conv = @benchmark convolve_step($q, $setup_p2p; params=$constants, ϵ=$ϵ)
    convolution[i, 1] = minimum(b_conv.times)

    # Line to point
    b_conv = @benchmark convolve_step($q, $setup_l2p; params=$constants, ϵ=$ϵ)
    convolution[i, 2] = minimum(b_conv.times)

    # Line to line
    b_conv = @benchmark convolve_step($q, $setup_l2l; params=$constants, ϵ=$ϵ)
    convolution[i, 3] = minimum(b_conv.times)
end