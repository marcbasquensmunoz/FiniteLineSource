using FiniteLineSource
using FiniteLineSource: PointSource
using BenchmarkTools
using SparseArrays
using FMM3D

ϵ = 1e-6
Δt = 3600.
# Nt = 8760*20
Nt = 50000

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

B = 5.

positions = [PointSource(B*i, B*j, B*k, rb) for i in 0:4 for j in 0:4 for k in 0:4 ]
q = 1e5 .* ones(length(positions),Nt)

# q = 1e5 .* repeat(vcat(ones(Nt - block.N[1]+1), zeros(block.N[1] -1)))

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = B)

# Block method
sources = targets = hcat([[p.x,p.y,p.z] for p in  positions]...)
res = zeros(size(targets)[2])

block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);

# FiniteLineSource.evolve_F!(q,1,8760*2, block)
FiniteLineSource.evolve_F!(q,1,Nt, block)
@time FiniteLineSource.fmm_evaluation!(res,sources,block; ffmeps = 1e-8)


# @time evolve!(Ib, q[:,1:8760*2], block)

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
    C  += convolve_step(q2, PointToPoint(r); params=constants)
end

# # Error
err = @. abs(Ib[1, :] - C)
##.

# #######################################
# # Performance analysis with non-history
# #######################################

# Inh = zeros(Nt)

# # Precomputation
# block = @btime prepare_containers(setup, positions, ϵ, Nt, constants);
# precomp = @btime precompute_parameters(setup, params=constants, ϵ=ϵ);

# # Simulation 
# @btime evolve!(Ib, q, block)
# @btime compute_integral_throught_history!(setup, I=Inh, q=q, precomp=precomp, params=constants)




# Rearrange block.Kranges so the ranges of each block is separated
Kb = push!([block.Kranges[j].start:block.Kranges[j+1].start-1 for j = 1:length(block.Kranges)-1], block.Kranges[end])

# Rearrange block.ζ so the ranges of each block is separated
ζb = [block.ζ[k] for k in Kb]

Ks1 = sparse((block.K_min .== 1))
Ks2 = sparse((block.K_min .== 2))
Ks3 = sparse((block.K_min .== 3))
Ks4 = sparse((block.K_min .== 4))
Ks5 = sparse((block.K_min .== 5))
Ks6 = sparse((block.K_min .== 6))
# Ks7 = sparse((block.K_min .== 7))


targets = hcat([[p.x,p.y,p.z] for p in  positions]...)
zk = complex(block.ζ[70])
charges = complex(ones(length(block.F[70,:])))
@time vals = hfmm3d(1e-12,zk,targets,charges=charges,targets = targets,pg=1)
imag.(vals.pot)