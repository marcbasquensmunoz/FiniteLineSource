using FiniteLineSource
using FiniteLineSource: PointSource
using BenchmarkTools

ϵ = 1e-6
Δt = 3600.
Nt = 8760

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

B = 5.
#positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
positions = [PointSource(B*i, B*j, B*k, rb) for i in 0:5 for j in 0:5 for k in 0:5]
q = ones(length(positions), Nt)

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = r)

# Block method
block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
Ib = zeros(length(positions), Nt)
@time evolve!(Ib, q, block)

C = zeros(Nt)
target = 1
# Convolution
for i in 1:length(positions)
    if target == i continue end
    r = FiniteLineSource.compute_distance_3D(positions[target], positions[i])
    C += convolve_step(q[i,:], PointToPoint(r = r); params=constants, ϵ=ϵ/length(positions))
end

# Error
err = @. abs(Ib[target, :] - C)
maximum(err)
