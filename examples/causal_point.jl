using SpecialFunctions
using FastGaussQuadrature
using LinearAlgebra
using Plots
using QuadGK
using FiniteLineSource

# Does not work for r̃ = 1
r = 0.7
α = 1e-6
kg = 3.
rb = 0.1

n = 10

Δt = 3600.
Nt = 1

h_ps(r, α, kg, t) = erfc(r/sqrt(4t*α)) / (4*π*r*kg)

#########################
Δt̃ = Δt*α/rb^2
r̃ = r/rb 

ϵ = 1e-12
t_ϵ = (r / sqrt(4α) / erfcinv(4*π*r*kg*ϵ))^2
Nt_prime = floor(t_ϵ/Δt)

#Q = [sin(2π /24 * i) + sin(2π / 8760 * i) for i in 1:Nt]
Q = [1 for i in 1:Nt]

res = 0.

if Nt_prime == 0
    I = zeros(length(Q))
    setup = PointToPoint(r=r)
    params = Constants(Δt = Δt, α=α, rb=rb, kg=kg)
    precomp = @time precompute_parameters(setup, params=params)
    @time compute_inFtegral_throught_history!(setup, I=I, q=Q, precomp=precomp, params=params)
    res = I[length(Q)]
else

    a = 0.
    b = sqrt(-log(ϵ/maximum(Q)) / (Nt_prime*Δt̃))

    guide(ζ) = exp(-ζ^2*Nt_prime*Δt̃) * sin(r̃*ζ)
    I, E, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)

    sort!(segbuf, by=x->x.a)
    n_seg = length(segbuf)

    N = n_seg*(n+1)

    ζ = zeros(N)
    W = zeros(N)
    F = zeros(N)
    x, w = gausslegendre(n+1)

    for (i, segment) in enumerate(segbuf)
        m = (segment.b-segment.a)/2
        c = (segment.b+segment.a)/2 
        @. ζ[(i-1)*(n+1)+1:i*(n+1)] = m*x + c
        @. W[(i-1)*(n+1)+1:i*(n+1)] = m * w
    end

    #########################
    # Non-history method
    #########################

    C = 1 / (2π^2*kg)
    expt = @. exp(-ζ^2*Δt̃)
    fact = @. (1 - expt) * sin(r̃*ζ) / (ζ * r)
    for q in Q
        @. F = expt * F + q * fact
    end

    ############################
    # Omit the loads that don't have time to reach the target
    ############################

    @. F = F * exp(-ζ^2*Nt_prime*Δt̃)

    res = C * dot(F, W)
end

error = res - h_ps(r, α, kg, Δt*(Nt+Nt_prime))
@show error
@show N

plot(ζ, F)

