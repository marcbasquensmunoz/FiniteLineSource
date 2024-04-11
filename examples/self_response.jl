using FiniteLineSource
using QuadGK
using SpecialFunctions
using Cubature

# Example of computation of self response
t = 3600.
rb = 0.1
H = 10.
α = 10^-6
kg = 3.

# Using method
@btime compute_self_response(3000, t=t, rb=rb, H=H, α=α, kg=kg)

# Integrating directly with representative
@btime begin 
params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=rb)
h_mean_sts, r_min, r_max = mean_sts_evaluation(params)
quadgk(r -> h_mean_sts(r) * point_step_response(t, r, α, kg), r_min, r_max)
end

# Integrating directly the double integral
@btime hcubature(z -> point_step_response(t, sqrt(rb^2 + (z[1]-z[2])^2), α, kg) / H, [0., 0.], [H, H])
@btime quadgk(z -> quadgk(zp -> point_step_response(t, sqrt(rb^2 + (z-zp)^2), α, kg) / H, 0., H)[1], 0., H)