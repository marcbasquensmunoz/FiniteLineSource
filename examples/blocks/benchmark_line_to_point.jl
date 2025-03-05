using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools
using Parameters

ϵ = 1e-6
Δt = 3600.
Nt = 10000

α = 1e-6
kg = 3.
rb = 0.1

Ds = 0.
Hs = 100.
Dt = 100.
Ht = 100.
B = 10.

q = [1. for t in 1:Nt]

bh_positions = [LineSource(x=0., y=0., D=Ds, H=Hs), LineSource(x=B, y=0., D=Dt, H=Ht)]
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb, line_points=[1, 1, 1, 1, 1] .* 500, line_limits=[0., 0.1, 0.3, 0.7, 0.9, 1.])

#####################################
# Error analysis
#####################################

# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(SegmentToPoint(D=0., H=0., z=0., σ=0.), bh_positions, ϵ, Nt, constants, containers);
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
setup = SegmentToPoint(D=Ds, H=Hs, z=Dt+Ht/2, σ=B)
C = convolve_step(q, setup; params=constants)

# Error
err = @. abs(Ib[1, :] - C)


#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(SegmentToPoint(D=0., H=0., z=0., σ=0.), bh_positions, ϵ, Nt, constants, containers);
precomp = @btime precompute_parameters(setup, params=constants);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)
