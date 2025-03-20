using FiniteLineSource
using BenchmarkTools

ϵ = 1e-6
Δt = 3600.
Nt = 1000

bn = 2
bm = 1
bl = 1

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

Δt̃ = Δt*α/rb^2

B = 0.7

q = [1. for t in 1:Nt]

positions = [(B*(i-1)^2, B*(j-1)^2, B*(k-1)^2) for i in 1:bn for j in 1:bm for k in 1:bl]


#####################################
# Error analysis
#####################################
setup = PointToPoint(r = B)

# Block method
block = prepare_containers(setup, positions, ϵ, Nt, constants, nothing);
Ib = zeros(length(positions), Nt)
evolve!(Ib, q, block)

# Convolution
C = convolve_step(q, setup; params=constants)

# Error
err = @. abs(Ib[1, :] - C)

#######################################
# Performance analysis with non-history
#######################################

Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(setup, positions, ϵ, Nt, constants, nothing);
precomp = @btime precompute_parameters(setup, params=params);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=params)
