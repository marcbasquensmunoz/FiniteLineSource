using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools
using Parameters

ϵ = 1e-6
Δt = 3600.
Nt = 20000

α = 1e-6
kg = 3.
rb = 0.1

Ds = 0.
Hs = 150.
Dt = 0.
Ht = 150.
B = 1.

q = [1. for t in 1:Nt]

bh_positions = [LineSource(x=0., y=0., D=Ds, H=Hs), LineSource(x=B, y=0., D=Dt, H=Ht)]
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

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
C = @time convolve_step(q, setup; params=constants)

# Error
err = @. abs(Ib[1, :] - C)


#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
precomp = @btime precompute_parameters(setup, params=constants);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)


Profile.clear()
@profile prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
pprof()

Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=0.5 prepare_containers(setup, bh_positions, ϵ, Nt, constants, containers);
PProf.Allocs.pprof(from_c=false)


line(h, N) = N <= 0 ? 0 : 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * abs(erf(s*(z-D-H/2+h)) - erf(s*(z-D-H/2-h))), 1/sqrt(4*α*Δt*N), Inf, atol=ϵ)[1]
NN = 1:2000
plot(NN, line.(20, NN))
Nrl = Int(floor(find_zero(N -> line(h, N) - ϵ, compute_N(d, ϵ, constants))))
