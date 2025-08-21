using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools
using Parameters

ϵ = 1e-8
Δt = 3600.
Nt = 8760*20

α = 1e-6
kg = 3.
rb = 0.1

Ds = 0.
Hs = 100.
Dt = 0.
Ht = 100.
B = 1.

q = hcat(
    [1. for t in 1:Nt], 
    [1. for t in 1:Nt], 
)'

bh_positions = [LineSource(x=0., y=0., D=Ds, H=Hs, rb=rb), LineSource(x=B, y=0., D=Dt, H=Ht, rb=rb)]
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

#####################################
# Error analysis
#####################################

setup = SegmentToPoint(D=Ds, H=Hs, z=Dt+Ht/2, σ=B, image_strength = -1.)
image_setup = SegmentToPoint(D=-Ds-Hs, H=Hs, z=Dt+Ht/2, σ=B)
sr_setup = SegmentToPoint(D=Ds, H=Hs, z=Ds+Hs/2, σ=rb)
sr_image_setup = SegmentToPoint(D=-Ds-Hs, H=Hs, z=Ds+Hs/2, σ=rb)

# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(setup, bh_positions, ϵ/length(bh_positions), Nt, constants, containers, Q=maximum(q), compute_self_response=true);
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
C_2to1 = @time convolve_step(q[2, :], setup; params=constants)
C_2to1_image = @time convolve_step(-q[2, :], image_setup; params=constants)
C_1sr = @time convolve_step(q[1, :], sr_setup; params=constants)
C_1sr_image = @time convolve_step(-q[1, :], sr_image_setup; params=constants)

C1 = C_2to1 + C_2to1_image + C_1sr + C_1sr_image

# Error
err = @. abs(Ib[1, :] - C1)

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
