using SpecialFunctions
using LinearAlgebra
using QuadGK
using DataStructures
using Parameters
using FastGaussQuadrature
using FiniteLineSource
using Roots

####### Compute the response of all points on all other points (excluding self-responses)

α = 1e-6
kg = 3.
rb = 0.1

n = 10
ϵ = 1e-12

Nt = 4
Δt = 3600.
Δt̃ = Δt*α/rb^2

params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

positions = [(0., 0., 0.), (1.1, 0., 0.), (0., 0.9, 0.), (1., 1., 0.), (10., 0., 0.)]
Q = [[1. for t in 1:Nt] for _ in positions]

compute_distance(x1, x2) = sqrt((x1[1]-x2[1])^2+(x1[2]-x2[2])^2+(x1[3]-x2[3])^2)
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
R = reduce(vcat, [[compute_distance(positions[i], positions[j]) for j in i+1:length(positions)] for i in 1:length(positions)])

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
    ζ[i], W[i] = compute_ζ_points(Ni, No, ϵ, maximum(maximum(Q)), n, params)
end

F = [[zeros(length(ζ[i])) for i in eachindex(N)] for _ in positions]

#####################################
# Precomputation
#####################################
expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(ζ)]
expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(ζ)]
expNout = [i == length(ζ) ? zeros(length(ζ[i])) : @. exp(-ζ[i]^2*(N[i+1])*Δt̃) for i in eachindex(ζ)]

load_delays = [[Queue{Float64}() for _ in eachindex(F[i])] for i in eachindex(F)]
load_buffers = [Queue{Float64}() for _ in eachindex(F)]

for j in 1:N[1]
    for buffer in load_buffers
        enqueue!(buffer, 0.)
    end
end
for source in load_delays
    for i in 1:length(source)-1
        for j in 1:N[i+1]-N[i]
            enqueue!(source[i], 0.)
        end
    end
end

#####################################
# Evolution in time
#####################################
for (i, source) in enumerate(Q)
    for q in source
        enqueue!(load_buffers[i], q)
        current_q = dequeue!(load_buffers[i])
        for j in eachindex(F[i])
            qin = current_q
            enqueue!(load_delays[i][j], qin)
            qout = dequeue!(load_delays[i][j])
            current_q = qout
            @. F[i][j] = expt[j] * F[i][j] + qin * expNin[j] - qout * expNout[j]
        end
    end
end

#####################################
# Reconstruction of integrand and integration
#####################################
C = 1 / (2π^2*kg)
In = zeros(length(positions))

for i in eachindex(positions)
    posR = compute_distance.(Ref(positions[i]), filter(x -> x != positions[i], positions))

    for r in posR
        r̃ = r/rb
        Nrpos = compute_N(r, ϵ, params)
        minN = findlast(n -> Nrpos >= n, N)
        for j in minN:length(F[i])
            In[i] += C * dot(F[i][j] .* (1 .- expt[j]) .* sin.(ζ[j] .* r̃) ./ (ζ[j] .* r), W[j])
        end
    end
end

#=
In = zeros(length(R))
for i in eachindex(R)
    r = R[i]
    r̃ = r/rb
    minN = findlast(n -> Nr[i] >= n, N)
    for j in minN:length(F)
        In[i] += C * dot(F[j] .* (1 .- expt[j]) .* sin.(ζ[j] .* r̃) ./ (ζ[j] .* r), W[j])
    end
end
=#

#####################################
# Analysis
#####################################
Ntot = sum(map(x->length(x), ζ[1:end]))


#####################################
# Validation
#####################################
using FiniteLineSource

I_real = zeros(length(positions))
for i in eachindex(positions)
    for j in eachindex(positions)
        i == j && continue
        setup = PointToPoint(r = compute_distance(positions[i], positions[j]))
        precomp = precompute_parameters(setup, params=params)
        Itemp = zeros(Nt)
        @views compute_integral_throught_history!(setup, I=Itemp, q=Q[i], precomp=precomp, params=params)
        I_real[i] += Itemp[end]
    end
end
err = @. In - I_real
