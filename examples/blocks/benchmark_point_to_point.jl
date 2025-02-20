using FiniteLineSource: Constants, precompute_parameters, compute_integral_throught_history!, PointToPoint, prepare_containers_ptp, evolve_ptp!
using BenchmarkTools

ϵ = 1e-6
Δt = 3600.
Nt = 100000

bn = 2
bm = 1
bl = 1

α = 1e-6
kg = 3.
rb = 0.1
params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

Δt̃ = Δt*α/rb^2

B = 5.

q = [1. for t in 1:Nt]

positions = [(B*(i-1)^2, B*(j-1)^2, B*(k-1)^2) for i in 1:bn for j in 1:bm for k in 1:bl]


#####################################
# Error analysis
#####################################

# Block method
N, block = prepare_containers_ptp(positions, ϵ, Nt, params);
Ib = zeros(length(positions), Nt)
evolve_ptp!(Ib, q, block)

# Original non-history
Inh = zeros(Nt)
setup = PointToPoint(r = B)
precomp = precompute_parameters(setup, params=params)
compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=params)

# Error
err = @. abs(Ib[1, :] - Inh)

#####################################
# Performance analysis
#####################################

# Precomputation
@btime prepare_containers_ptp(positions, ϵ, Nt, params);
@btime precompute_parameters(setup, params=params)

# Simulation 
@btime evolve_ptp!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=params)
