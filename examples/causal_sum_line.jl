using SpecialFunctions
using LinearAlgebra
using QuadGK
using DataStructures
using Parameters
using FastGaussQuadrature
using FiniteLineSource
using Roots

include("block_method.jl")


α = 1e-6
kg = 3.
rb = 0.1

D = 0.
H = 100.

n = 10
ϵ = 1e-12

Nl = 400
Nt = 100000
Δt = 3600.
Δt̃ = Δt*α/rb^2


params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [1. for t in 1:Nt]

# Evaluation points 
# DO NOT USE FOR SELF-RESPONSE
R = [1.]#[1., 2., 5., 10., 11., 13., 20., 30., 50., 100.]

#####################################
# Initialization
#####################################
#N = compute_N.(R, ϵ, Ref(params))
#Nmin = compute_N_line(minimum(R), Ref(D), Ref(H), Ref(D+H/2), Ref(ϵ), Ref(params))
#Nmax = compute_N_line(maximum(R), Ref(D), Ref(H), Ref(D+H/2), Ref(ϵ), Ref(params))

Nmin = compute_N(minimum(R), ϵ, params)
Nmax = compute_N(sqrt(maximum(R)^2+(H/2)^2), ϵ, params)

ratio = 5
Ncurrent = Nmax / ratio
N = zeros(Int, 0)

while Ncurrent >= Nmin
    push!(N, Int(floor(Ncurrent)))
    Ncurrent /= ratio
end
!(Nmin in N) && push!(N, Nmin)
sort!(N)
Nr = compute_N_line.(R, Ref(D), Ref(H), Ref(D+H/2), Ref(ϵ), Ref(params))

for (na, nb) in zip(N[1:end-1], N[2:end])
    if isempty(filter(n -> n in na:nb, Nr)) 
        deleteat!(N, findfirst(x->x==na, N))
    end
end

ζ = [zeros(0) for _ in eachindex(N)]
W = [zeros(0) for _ in eachindex(N)]
for i in eachindex(N)
    Ni = N[i]
    No = i == length(N) ? 3*N[i] : N[i+1]
    ζ[i], W[i] = compute_ζ_points_line(Ni, No, ϵ, maximum(q), n, D, H, D+H/2, params)
end

F = [zeros(length(ζ[i])) for i in eachindex(N)]

#####################################
# Precomputation
#####################################
expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(F)]
expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(F)]
# Does this go with a N[i+1]+1 ??? 
expNout = [i == length(F) ? zeros(length(F[i])) : @. exp(-ζ[i]^2*(N[i+1])*Δt̃) for i in eachindex(F)]

load_delays = [Queue{Float64}() for _ in eachindex(F)]
load_buffer = Queue{Float64}()

for j in 1:N[1]
    enqueue!(load_buffer, 0.)
end
for i in 1:length(load_delays)-1
    for j in 1:N[i+1]-N[i]
        enqueue!(load_delays[i], 0.)
    end
end

#####################################
# Evolution in time
#####################################
for qt in q
    enqueue!(load_buffer, qt)
    current_q = dequeue!(load_buffer)

    for i in eachindex(F)
        qin = current_q
        enqueue!(load_delays[i], qin)
        qout = dequeue!(load_delays[i])
        current_q = qout
        @. F[i] = expt[i] * F[i] + qin * expNin[i] - qout * expNout[i]
    end
end

#####################################
# Analysis
#####################################
Ntot = sum(map(x->length(x), ζ[1:end]))


#####################################
# Complete F
#####################################
C = 1 / (2π^2*kg)
for i in eachindex(F)
    @. F[i] = F[i] * (1 - expt[i]) / ζ[i]
end

#####################################
# Line integral with HMatrices
#####################################
using HMatrices
using StaticArrays
const Point3D = SVector{3,Float64}

z, wz = gausslegendre(Nl)
Z = @. H/2 * z + D + H/2
WZ = @. H/2 * wz

sources = [Point3D(0., 0., zz) for zz in Z]
targets = [Point3D(r, 0., D+H/2) for r in R]

function kernel(k, rb)
    function G(x,y)
        d = norm(x-y)
        sin(k*d/rb) / d
    end
    return G
end

HM = [[assemble_hmatrix(KernelMatrix(kernel(ζζ, rb), sources, targets); atol=ϵ) for ζζ in ζ[i]] for i in eachindex(ζ)]

#=
R_block_index = [findlast(x -> x <= n, N) for n in Nr]
Ipart = zeros(length(R))
for j in length(F):-1:1
    update = 1:findlast(x->x==j, R_block_index)
    for k in eachindex(ζ[j])
        @views Ipart[update] .+= ((C * F[j][k] * W[j][k]) .* (WZ' * HM[j][k])')[update]
    end
end
=#

Istp = zeros(length(R))

function blocks_for_line(i, z_eval, Z, eval_N)
    compute_N_point_in_line(z) = compute_N(sqrt(R[i]^2+(z_eval-z)^2), ϵ, params)
    Ns = compute_N_point_in_line.(Z)
    map(n -> findlast(m -> m <= n, eval_N), Ns)
end

Nm = [blocks_for_line(i, D+H/2, Z, N) for i in 1:length(R)]
Nmtot = [minimum(reshape(reduce(vcat,Nm), (length(Z), length(R)))[i,:]) for i in 1:length(Z)]
R_block_index = [findlast(x -> x <= n, N) for n in Nr]
for j in length(F):-1:1
    line_points = findall(x-> x <= j, Nmtot)
    AA = findlast(x->x==j, R_block_index)
    update = 1:(isnothing(AA) ? length(R) : AA)
    WWZ = zeros(length(WZ))
    WWZ[line_points] = WZ[line_points]
    for k in eachindex(ζ[j])
        @views Istp[update] .+= ((C * F[j][k] * W[j][k]) .* (WWZ' * HM[j][k])')[update]
    end
end

#####################################
# Validation
#####################################
using FiniteLineSource

STP_test = zeros(length(R))
for i in eachindex(R)
    σ = R[i]
    STP_test[i] = quadgk(zp -> erfc(sqrt(σ^2 + (zp - D-H/2)^2) / sqrt(4α * Δt * Nt))/(4*π*kg*sqrt(σ^2 + (zp -  D-H/2)^2)), D, D+H)[1]
end

err = @. abs(Istp - STP_test)
