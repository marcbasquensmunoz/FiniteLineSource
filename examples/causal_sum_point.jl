using SpecialFunctions
using LinearAlgebra
using QuadGK
using DataStructures
using Parameters
using FastGaussQuadrature
using FiniteLineSource
using Roots

α = 1e-6
kg = 3.
rb = 0.1

n = 10
ϵ = 1e-12

Nt = 100
Δt = 3600.
Δt̃ = Δt*α/rb^2


params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [1. for t in 1:Nt]

# Evaluation points 
# DO NOT USE FOR SELF-RESPONSE
R = [1.]


"""
Compute the number of steps N that can be skipped for a given distance r
"""
function compute_N(r, ϵ, params::Constants) 
    @unpack Δt, α, kg = params
    Int(floor((r / sqrt(4α) / erfcinv(4*π*r*kg*ϵ))^2/Δt))
end
"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points(N, No, ϵ, Q, n, params::Constants)
    @unpack Δt, α, rb = params
    Δt̃ = Δt*α/rb^2

    heatwave(r) = erfc(r/sqrt(4α*N*Δt)) / r - ϵ
    r = find_zero(heatwave, sqrt(N))
    a = 0.
    b = sqrt(-log(ϵ/Q) / (N*Δt̃))
    r̃ = r/rb

    guide(ζ) = (exp(-ζ^2*N*Δt̃) + 100 * exp(-ζ^2*No*Δt̃)) * sin(r̃*ζ) / (r*ζ) * (1 - exp(-ζ^2*Δt̃))
    #guide(ζ) = exp(-ζ^2*N*Δt̃) * sin(r̃*ζ) / (r*ζ) * (1 -  exp(-ζ^2*Δt̃))
    _, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    sort!(segbuf, by=x->x.a)
    n_seg = length(segbuf)

    Nζ = n_seg*n

    ζ = zeros(Nζ)
    W = zeros(Nζ)
    x, w = gausslegendre(n)

    for (i, segment) in enumerate(segbuf)
        m = (segment.b-segment.a)/2
        c = (segment.b+segment.a)/2 
        @. ζ[(i-1)*n+1:i*n] = m*x + c
        @. W[(i-1)*n+1:i*n] = m * w
    end

    return ζ, W
end

#####################################
# Initialization
#####################################
#N = compute_N.(R, ϵ, Ref(params))
Nmin = compute_N(minimum(R), ϵ, params)
Nmax = compute_N(maximum(R), ϵ, params)

ratio = 5
Ncurrent = Nmax / ratio
N = zeros(Int, 0)

while Ncurrent >= Nmin
    push!(N, Int(floor(Ncurrent)))
    Ncurrent /= ratio
end
!(Nmin in N) && push!(N, Nmin)
sort!(N)
Nr = compute_N.(R, ϵ, Ref(params))

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
    ζ[i], W[i] = compute_ζ_points(Ni, No, ϵ, maximum(q), n, params)
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
# Reconstruction of integrand and integration
#####################################
C = 1 / (2π^2*kg)
In = zeros(length(R))
for i in eachindex(R)
    r = R[i]
    r̃ = r/rb
    minN = findlast(n -> Nr[i] >= n, N)
    for j in minN:length(F)
        In[i] += C * dot(F[j] .* (1 .- expt[j]) .* sin.(ζ[j] .* r̃) ./ (ζ[j] .* r), W[j])
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

I_test = I_real[Nt, :]

err = @. In - I_test
