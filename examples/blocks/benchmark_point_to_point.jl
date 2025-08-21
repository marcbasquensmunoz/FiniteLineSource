using FiniteLineSource
using FiniteLineSource: PointSource

ϵ = 1e-6
Δt = 3600.
Nt = 8760*20

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

B = 1.

q = hcat(
    [1. for t in 1:Nt], 
    [1. for t in 1:Nt], 
)'

positions = [PointSource(0., 0., 0., rb), PointSource(B, 0., 0., rb)]

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = B)

# Block method
block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
Ib = zeros(length(positions), Nt)
evolve!(Ib, q, block)

# Convolution
C_int = convolve_step(q[2,:], setup; params=constants)
C_sr = convolve_step(q[1,:], PointToPoint(r=rb); params=constants)

C = C_int #+ C_sr
# Error
err = @. abs(Ib[1, :] - C)

#######################################
# Performance analysis with non-history
#######################################

Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(setup, positions, ϵ, Nt, constants);
precomp = @btime precompute_parameters(setup, params=constants, ϵ=ϵ);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)
