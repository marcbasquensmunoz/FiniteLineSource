using SpecialFunctions
using LinearAlgebra
using QuadGK
using DataStructures
using Parameters
using FastGaussQuadrature
using FiniteLineSource
using FiniteLineSource: compute_N, compute_ζ_points
using Roots

α = 1e-6
kg = 3.
rb = 0.1

n = 10
ϵ = 1e-6

Nt = 1000
Δt = 3600.
Δt̃ = Δt*α/rb^2


params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [1. for t in 1:Nt]

# Evaluation points 
# DO NOT USE FOR SELF-RESPONSE
R = [5.]


#####################################
# Initialization
#####################################
Nr = [compute_N(r, ϵ, params) for r in R]
N = choose_blocks(Nr, Nt)

ζ = [zeros(0) for _ in eachindex(N)]
W = [zeros(0) for _ in eachindex(N)]
for i in eachindex(N)
    Ni = N[i]
    No = i == length(N) ? 3*N[i] : N[i+1]
    ζ[i], W[i] = compute_ζ_points(Ni, No, ϵ, maximum(q), n, params)
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
In = zeros(Nt, length(R))
C = 1 / (2π^2*kg)

for (t, qt) in enumerate(q)
    enqueue!(load_buffer, qt)
    current_q = dequeue!(load_buffer)

    for i in eachindex(F)
        qin = current_q
        enqueue!(load_delays[i], qin)
        qout = dequeue!(load_delays[i])
        current_q = qout
        @. F[i] = expt[i] * F[i] + qin * expNin[i] - qout * expNout[i]
    end

    for i in eachindex(R)
        r = R[i]
        r̃ = r/rb
        minN = findlast(n -> Nr[i] >= n, N)
        for j in minN:length(F)
            In[t, i] += C * dot(F[j] .* (1 .- expt[j]) .* sin.(ζ[j] .* r̃) ./ (ζ[j] .* r), W[j])
        end
    end
end

#####################################
# Analysis
#####################################
Ntot = sum(map(x->length(x), ζ[1:end]))


#####################################
# Validation
#####################################
using FiniteLineSource

I_real = zeros(length(q), length(R))

for i in eachindex(R)
    setup = PointToPoint(r = R[i])
    precomp = precompute_parameters(setup, params=params)
    @views compute_integral_throught_history!(setup, I=I_real[:, i], q=q, precomp=precomp, params=params)
end


err = @. In - I_real
