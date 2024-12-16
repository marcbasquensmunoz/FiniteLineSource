using SpecialFunctions
using LinearAlgebra
using QuadGK
using DataStructures
using Parameters
using FastGaussQuadrature
using FiniteLineSource
using Roots
using Bessels

include("block_method.jl")

α = 1e-6
kg = 3.
rb = 0.1

D = 0.
H = 100.

n = 10
Nl = 500
ϵ = 1e-8

Nt = 1000
Δt = 3600.
Δt̃ = Δt*α/rb^2

z_eval = 10.


params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [1. for t in 1:Nt]

# Evaluation points 
# DO NOT USE FOR SELF-RESPONSE
R = [20.]#[1., 3., 5., 10.]#[1., 2., 5., 10., 11., 13., 20., 30., 50., 100.]

#####################################
# Initialization
#####################################
Nr = compute_N_line.(R, Ref(D), Ref(H), Ref(z_eval), Ref(ϵ), Ref(params))
#Nr = compute_N.(R, Ref(ϵ), Ref(params))

#N = [3, 300, 1000]#
N = choose_blocks(Nr, p = 10)

ζ = [zeros(0) for _ in eachindex(N)]
W = [zeros(0) for _ in eachindex(N)]

for i in eachindex(N)
    Ni = N[i]
    No = i == length(N) ? 3*N[i] : N[i+1]
    ζ[i], W[i] = compute_ζ_points_line(Ni, No, ϵ, maximum(q), n, D, H, z_eval, params)
end

F = [zeros(length(ζ[i])) for i in eachindex(N)]

#####################################
# Precomputation
#####################################
expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(F)]
expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(F)]
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
# Line integral 
#####################################
using HMatrices
using StaticArrays
const Point3D = SVector{3,Float64}


sources = [Point3D(0., 0., 0.)]
targets = [Point3D(r, 0., z_eval) for r in R]

xb, wb = gausslegendre(Nl+1) 
Pb = [wb[s] * Pl(xb[s], k) for k = 0:Nl, s = 1:Nl+1]
Xb = zeros(Nl+1)
aux = zeros(Nl+1)

function G(x, y, k)
    setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2))
    compute_kernel_line(k, rb, setup; x=xb, X=Xb, P=Pb, f=aux, atol=ϵ)   
end

HM = @time [[[G(source, target, ζζ) for source in sources, target in targets] for ζζ in ζ[i]] for i in eachindex(ζ)]

R_block_index = [findlast(x -> x <= n, N) for n in Nr]
Istp = zeros(length(R))
for j in length(F):-1:1
    update = findall(x -> x <= j, R_block_index)
    upM = diagm(ones(length(R)))
    for k in eachindex(ζ[j])
        @views Istp[update] .+= ((C * F[j][k] * W[j][k]) * HM[j][k] * upM)'[update]
    end
end


#####################################
# Validation
#####################################

STP_test = zeros(length(R))
for i in eachindex(R)
    σ = R[i]
    STP_test[i] = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z_eval)^2) / sqrt(4α * Δt * Nt))/(4*π*kg*sqrt(σ^2 + (zp - z_eval)^2)), D, D+H, atol = ϵ)[1]
end

err = @. abs(Istp - STP_test)
