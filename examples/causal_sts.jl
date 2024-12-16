using SpecialFunctions
using FastGaussQuadrature
using LinearAlgebra
using Plots
using QuadGK
using FiniteLineSource
using FiniteLineSource: h_mean_lims, h_mean_sts
    
α = 1e-6
kg = 3.
rb = 0.1

D1 = 0.
D2 = 0.
H1 = 100.
H2 = 100.
σ = 1.

n = 6
nr = 6

Δt = 3600.
Nt = 100#8760

ierf(x) = x*erf(x) - (1 - exp(-x^2)) / sqrt(π)
function h_sts(D1, H1, D2, H2, σ, α, kg, t) 
    I(s) = ierf((D2 - D1 + H2)*s) + ierf((D2 - D1 - H1)*s) - ierf((D2 - D1)*s) - ierf((D2 - D1 + H2 - H1)*s) 
    return 1/(4π*kg) * quadgk(s -> exp(-σ^2 * s^2) / s^2 * I(s), 1/sqrt(4α*t), Inf, atol=eps())[1]
end

sts_params = MeanSegToSegEvParams(SegmentToSegment(σ=σ, D1=D1, D2=D2, H1=H1, H2=H2))
r_min, r_max = h_mean_lims(sts_params) 
h_sts(r) = h_mean_sts(r, sts_params)

#########################
Δt̃ = Δt*α/rb^2 

ϵ = 1e-12
t_ϵ = (r_min / sqrt(4α) / erfcinv(4*π*σ*kg*ϵ))^2
Nt_prime = floor(t_ϵ/Δt)

Q = [1 for i in 1:Nt]

a = 0.
b = sqrt(-log(ϵ/maximum(Q)) / (Nt_prime*Δt̃))

guide(ζ) = exp(-ζ^2*Nt_prime*Δt̃) * sin(r_max/rb*ζ)
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

guide_r(r) = erf(r/rb) / r *rb * h_sts(r)
I, E, segbuf_r = quadgk_segbuf(guide_r, r_min, r_max, atol=ϵ)

sort!(segbuf_r, by=x->x.a)
nr_seg = length(segbuf_r)
Nr = nr_seg*(n+1)

R = zeros(Nr)
WR = zeros(Nr)

for (i, segment) in enumerate(segbuf_r)
    m = (segment.b-segment.a)/2
    c = (segment.b+segment.a)/2 
    @. R[(i-1)*(n+1)+1:i*(n+1)] = m * x + c
    @. WR[(i-1)*(n+1)+1:i*(n+1)] = m * w
end


#########################
# Non-history method
#########################

F = zeros(N, Nr)
GG = zeros(N, Nr)

expt = @. exp(-ζ^2*Δt̃)

for i in 1:Nr
    r = R[i]
    r̃ = r / rb
    @views @. GG[:, i] = (1 - expt) * sin(r̃*ζ) / (ζ * r)
end

C = 1 / (2π^2*kg)
for q in Q
    for i in 1:Nr
        @views @. F[:, i] = expt * F[:, i] + q * GG[:, i]
    end
end

############################
# Omit the loads that don't have time to reach the target
############################

for i in 1:Nr
    @views @. F[:, i] = F[:, i] * exp(-ζ^2*Nt_prime*Δt̃)
end

error = C * (W' * F * WR) - h_sts(D1, H1, D2, H2, σ, α, kg, Δt*(Nt+Nt_prime))
@show error
@show N, Nr
