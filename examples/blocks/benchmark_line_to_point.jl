using FiniteLineSource
using BenchmarkTools
using Parameters

ϵ = 1e-6
Δt = 3600.
Nt = 1000

bn = 1
bm = 2

α = 1e-6
kg = 3.
rb = 0.1

D = 0.
H = 100.
B = 5.

q = [1. for t in 1:Nt]

@with_kw struct LineSource{T <: Number} @deftype T
    x
    y
    D
    H
end

bh_positions = [LineSource(x=B*(i-1)^2, y=B*(j-1)^2, D=D, H=H) for i in 1:bn for j in 1:bm]
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb, line_points=[1, 1, 1, 1, 1] .* 500, line_limits=[0., 0.1, 0.3, 0.7, 0.9, 1.])

#####################################
# Error analysis
#####################################

# Block method
nmodel = FiniteLineSource.load_nmodel()
block = @time prepare_containers_ltp(bh_positions, ϵ, Nt, constants, nmodel);
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
setup = SegmentToPoint(D=D, H=H, z=D+H/2, σ=B)
C = convolve_step(q, setup; params=constants)

# Error
err = @. abs(Ib[1, :] - C)


#######################################
# Performance analysis with non-history
#######################################
Inh = zeros(Nt)

# Precomputation
block = @btime prepare_containers_ltp(bh_positions, ϵ, Nt, constants, nmodel);
precomp = @btime precompute_parameters(setup, params=constants);

# Simulation 
@btime evolve!(Ib, q, block)
@btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)
