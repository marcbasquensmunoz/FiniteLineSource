using FiniteLineSource
using FiniteLineSource: PointSource
using BenchmarkTools

ϵ = 1e-6
Δt = 3600.
Nt = 8760*20

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

r = 1.

positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
q = ones(length(positions), Nt)

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = r)

# Block method
block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
Ib = zeros(length(positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
C = convolve_step(q[2,:], setup; params=constants)

# Error
err = @. abs(Ib[1, :] - C)
maximum(err)
