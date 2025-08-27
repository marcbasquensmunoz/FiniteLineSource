using FiniteLineSource
using FiniteLineSource: PointSource
using BenchmarkTools
using SparseArrays
using FMM3D

ϵ = 1e-6
Δt = 3600.
Nt = 50000

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

B = 5.

positions = [PointSource(B*i, B*j, B*k, rb) for i in 0:2 for j in 0:2 for k in 0:2]
q = ones(length(positions), Nt)

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = B)

# Block method
sources = targets = hcat([[p.x,p.y,p.z] for p in  positions]...)
res = zeros(size(targets)[2])

block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
FiniteLineSource.evolve_F!(q, 1, Nt, block)
@time FiniteLineSource.fmm_evaluation!(res, sources, block; fmmeps = 1e-8)

block2 = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
Ib = zeros(length(positions), Nt)
@time evolve!(Ib, q[:,1:Nt], block2)

abs.(res .- Ib[:,Nt])
maximum(abs.(res .- Ib[:,Nt]))
##.

# # Convolution

q2 = vcat(ones(Nt - block.N[1]+1), zeros(block.N[1] -1))
C  = zeros(Nt) 
for pos in positions[2:end]
    r  = FiniteLineSource.compute_distance_3D(positions[1],pos)
    C  += convolve_step(q2, PointToPoint(r=r); params=constants)
end

# # Error
err = @. abs(Ib[1, :] - C)
maximum(err)

res[1] - C[end]
##.


# Rearrange block.ζ so the ranges of each block is separated
ζb = [block.ζ[k] for k in block.ranges]

Ks1 = sparse((block.K_min .== 1))
Ks2 = sparse((block.K_min .== 2))
Ks3 = sparse((block.K_min .== 3))
Ks4 = sparse((block.K_min .== 4))
Ks5 = sparse((block.K_min .== 5))
Ks6 = sparse((block.K_min .== 6))
