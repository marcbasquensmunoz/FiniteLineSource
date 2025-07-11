using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools


ϵ = 1e-4
Nt = 8760
B = 1.

D1 = 0.
D2 = 50.
H1 = 100.
H2 = 100.

α = 1e-6
kg = 3.
rb = 0.0575
Δt = 3600.
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [.5 for t in 1:Nt] 
#q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:Nt]

bh_positions = [LineSource(x=0., y=0., D=D1, H=H1, rb=rb), LineSource(x=B, y=0., D=D2, H=H2, rb=rb)]

#####################################
# Error analysis
#####################################

setup = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=B)
# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);#, Q=maximum(q));
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, ones(2) * q', block)

# Convolution
C = convolve_step(q, SegmentToSegmentOld(setup); params=constants, ϵ=ϵ/100)
C_sr = convolve_step(q, SegmentToSegmentOld(D1=D1, H1=H1, D2=D2, H2=H2, σ=rb); params=constants, ϵ=ϵ/100)
Ct = C+C_sr
# Error
err = @. abs(Ib[1, :] - Ct)

maximum(err)

#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
precomp = @btime precompute_parameters(SegmentToSegmentOld(setup), params=constants, ϵ=ϵ);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)



setup = SegmentToSegment(D1=0., H1=100., D2=100., H2=100., σ=1.)

f(h) = FiniteLineSource.compute_N_for_line_to_line_range(h, setup, constants, ϵ, Nt, 134, 1.)
hh = 0.1:0.1:50
plot(hh, f.(hh))