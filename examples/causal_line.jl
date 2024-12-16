using SpecialFunctions
using FastGaussQuadrature
using LinearAlgebra
using Plots
using QuadGK
    
α = 1e-6
kg = 3.
rb = 0.1

D = 0.
H = 100.
z = 50.
σ = 1.

n = 6
nr = 6

Δt = 3600.
Nt = 1#8760

h_ps(r, α, kg, t) = erfc(r/sqrt(4t*α)) / (4*π*r*kg)
h_stp(D, H, σ, z, α, kg, t) = quadgk(ζ -> h_ps(sqrt(σ^2 + (z-ζ)^2), α, kg, t), D, D+H)[1]

#########################
Δt̃ = Δt*α/rb^2
r = sqrt(σ^2 + (H-z)^2)
r̃ = r/rb 

ϵ = 1e-12
t_ϵ = (σ / sqrt(4α) / erfcinv(4*π*σ*kg*ϵ))^2
Nt_prime = floor(t_ϵ/Δt)

Q = [1 for i in 1:Nt]

a = 0.
b = sqrt(-log(ϵ/maximum(Q)) / (Nt_prime*Δt̃))

guide(ζ) = exp(-ζ^2*Nt_prime*Δt̃) * sin(r̃*ζ)
I, E, segbuf = quadgk_segbuf(guide, a, b, atol=ϵ)

sort!(segbuf, by=x->x.a)
n_seg = length(segbuf)

N = n_seg*(n+1)

ζ = zeros(N)
W = zeros(N)
x, w = gausslegendre(n+1)

for (i, segment) in enumerate(segbuf)
    m = (segment.b-segment.a)/2
    c = (segment.b+segment.a)/2 
    @. ζ[(i-1)*(n+1)+1:i*(n+1)] = m*x + c
    @. W[(i-1)*(n+1)+1:i*(n+1)] = m * w
end

guide_r(zp) = erf(sqrt(σ^2 + (z-zp)^2)) / sqrt(σ^2 + (z-zp)^2)
I, E, segbuf_r = quadgk_segbuf(guide_r, D, D+H, atol=ϵ)

sort!(segbuf_r, by=x->x.a)
nr_seg = length(segbuf_r)
Nr = nr_seg*(n+1)

Z = zeros(Nr)
WZ = zeros(Nr)

for (i, segment) in enumerate(segbuf_r)
    m = (segment.b-segment.a)/2
    c = (segment.b+segment.a)/2 
    @. Z[(i-1)*(n+1)+1:i*(n+1)] = m * x + c
    @. WZ[(i-1)*(n+1)+1:i*(n+1)] = m * w
end


#########################
# Non-history method
#########################

@time begin
F = zeros(N, Nr)
GG = zeros(N, Nr)

expt = @. exp(-ζ^2*Δt̃)

for i in 1:Nr
    r = sqrt(σ^2 + (z - Z[i])^2)
    r̃ = r / rb
    @views @. GG[:, i] = (1 - expt) * sin(r̃*ζ) / (ζ * r)
end

C = 1 / (2π^2*kg)
for q in Q
    for i in 1:Nr
        @views @. F[:, i] = expt * F[:, i] + q * GG[:, i]
    end
end
end

############################
# Omit the loads that don't have time to reach the target
############################

for i in 1:Nr
    @views @. F[:, i] = F[:, i] * exp(-ζ^2*Nt_prime*Δt̃)
end

C * (W' * F * WZ)

error = C * (W' * F * WZ) - h_stp(D, H, σ, z, α, kg, Δt*(Nt+Nt_prime))
@show error
@show N, Nr
