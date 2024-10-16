using FiniteLineSource
using FiniteLineSource: adaptive_gk, h_mean_lims, h_mean_sts, adaptive_nodes_and_weights
using QuadGK
using SpecialFunctions
using Cubature
using LinearAlgebra


function test_classical(;t, rb, H, α, kg)
    I_fls(s) = 2*ierf(H*s)
    ierf(x) = x*erf(x) - 1/sqrt(π)*(1- exp(-x^2))
    quadgk(s -> exp(-rb^2*s^2) / s^2 * I_fls(s) / (4H*π*kg), 1/sqrt(4*α*t), Inf)
end

function test_respresentative(;t, rb, H, α, kg)
    params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=rb)
    r_min, r_max = h_mean_lims(params)
    quadgk(r -> h_mean_sts(r, params) * point_step_response(t, r, α, kg), r_min, r_max)
end

function adaptive(;t, rb, H, α, kg)
    params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=rb)
    r_min, r_max = h_mean_lims(params)
    f(r) = h_mean_sts(r, params) * point_step_response(t, r, α, kg)
    x, w = adaptive_gk(f, r_min, r_max)
    dot(f.(x), w)
end

function adaptive_2(;t, rb, H, α, kg)
    params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=rb)
    r_min, r_max = h_mean_lims(params)
    f(r) = h_mean_sts(r, params) * point_step_response(t, r, α, kg)
    x, w = adaptive_nodes_and_weights(f, r_min, r_max)
    dot(f.(x), w)
end

# Example of computation of self response
t = 3600.
rb = 0.1
H = 100.
α = 10^-6
kg = 3.

# Using method
@time compute_self_response(300, t=t, rb=rb, H=H, α=α, kg=kg)

# Integrating directly with representative h
@time test_respresentative(t=t, rb=rb, H=H, α=α, kg=kg)

# Integrating directly the double integral
#@btime hcubature(z -> point_step_response(t, sqrt(rb^2 + (z[1]-z[2])^2), α, kg) / H, [0., 0.], [H, H])
@time quadgk(z -> quadgk(zp -> point_step_response(t, sqrt(rb^2 + (z-zp)^2), α, kg) / H, 0., H)[1], 0., H)

# Or with the classical formula
@time test_classical(t=t, rb=rb, H=H, α=α, kg=kg)

# Adaptive discretization
@time adaptive(t=t, rb=rb, H=H, α=α, kg=kg)
@time adaptive_2(t=t, rb=rb, H=H, α=α, kg=kg)


#### THIS DOESNT WORK WITH BIG H