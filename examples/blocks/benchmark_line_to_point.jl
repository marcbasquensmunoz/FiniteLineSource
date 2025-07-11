using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools
using Parameters

ϵ = 1e-4
Δt = 3600.
Nt = 8760

α = 1e-6
kg = 3.
rb = 0.1

Ds = 0.
Hs = 150.
Dt = 0.
Ht = 150.
B = 1.46779926762207

q = [1. for t in 1:Nt]
q2 = [20*sin(2π*t/8760) + 5*sin(2π*t/24) + 5. for t in 1:Nt]

bh_positions = [LineSource(x=0., y=0., D=Ds, H=Hs, rb=rb), LineSource(x=B, y=0., D=Dt, H=Ht, rb=rb)]
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

#####################################
# Error analysis
#####################################

setup = SegmentToPoint(D=Ds, H=Hs, z=Dt+Ht/2, σ=B)

# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers, Q=maximum(q2));
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, ones(2) * q', block)
@time evolve!(Ib, hcat(q, q2)', block)

# Convolution
C = @time convolve_step(q, setup; params=constants)
C2 = @time convolve_step(q2, setup; params=constants)
C_sr = @time convolve_step(q, SegmentToPoint(D=Ds, H=Hs, z=Ds+Hs/2, σ=rb); params=constants)
C_sr_2 = @time convolve_step(q2, SegmentToPoint(D=Ds, H=Hs, z=Ds+Hs/2, σ=rb); params=constants)

Ct = C + C_sr
# Error
err = @. abs(Ib[1, :] - C_sr - C2)
err = @. abs(Ib[2, :] - C_sr_2 - C)


maximum(err)


#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
precomp = @btime precompute_parameters(setup, params=constants, ϵ=ϵ);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)


#######################################
# Non-history error analysis
#######################################
Inh = zeros(Nt)
precomp = precompute_parameters(setup, params=constants);
compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)
maximum(abs.(Inh-C))



FiniteLineSource.compute_N_for_line_range(0., setup, constants, ϵ, Nt, 10000, 1.)
