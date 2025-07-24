using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools


ϵ = 1e-6
Nt = 8760*20
B = 1.

D1 = 0.
D2 = 0.
H1 = 100.
H2 = 100.

α = 1e-6
kg = 3.
rb = 0.0575
Δt = 3600.
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = hcat(
    [1. for t in 1:Nt],
    [1. for t in 1:Nt]
)'

bh_positions = [LineSource(x=0., y=0., D=D1, H=H1, rb=rb), LineSource(x=B, y=0., D=D2, H=H2, rb=rb)]

#####################################
# Error analysis
#####################################

setup = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=B)
# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(setup, bh_positions, ϵ/length(bh_positions), Nt, constants, containers, Q=maximum(q));
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
C_int = @time convolve_step(q[2,:], SegmentToSegmentOld(setup); params=constants, ϵ=ϵ/100)
C_sr = @time convolve_step(q[1,:], SegmentToSegmentOld(D1=D1, H1=H1, D2=D2, H2=H2, σ=rb); params=constants, ϵ=ϵ/100)
C = C_int + C_sr
# Error
err = @. abs(Ib[1, :] - C)

maximum(err)


#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb, line_points=[50, 50], line_limits=[0., 0.5, 1.])

# Precomputation
block = @btime prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
precomp = @btime precompute_parameters(SegmentToSegmentOld(setup), params=constants, ϵ=ϵ);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)
