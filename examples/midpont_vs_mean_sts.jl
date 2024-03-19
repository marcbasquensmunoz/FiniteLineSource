using FiniteLineSource
using FiniteLineSource: point_step_response, segment_step_response, segment_mean_step_response
using Plots 
using QuadGK

rb = 0.1
σ = 5.
α = 10^-6
kg = 3.
t = 3600*24*30*12.

## Source line
D = 0.
H = 50.

## Target line
Dt = 20.
Ht = 20.

####

params = LevelSetParams(D1=D, H1=H, D2=Dt, H2=Ht, σ=σ, rb=rb)
h_mean_sts, r_mean = mean_sts_evaluation(params)

## Comparison of level set method with actual double integral
pred_mean = quadgk(r -> h_mean_sts(r)*point_step_response(t, rb*r, α, kg), r_mean[1]+10^-14, r_mean[end])[1]
exact_mean = segment_mean_step_response(t, σ, params.D1, params.H1, params.D2, params.H2, α, kg)
pred_mean - exact_mean

Plots.plot(r_mean, h_mean_sts.(r_mean))
####

params_mid = MidPointParams(D=D, H=H, z=Dt+Ht/2, σ=σ, rb=rb)
h_midpoint, r_mid = midpoint_evaluation(params_mid)

## Comparison of midpoint method with actual double integral
pred_mid = quadgk(r -> h_midpoint(r)* point_step_response(t, r*rb, α, kg), r_mid[1], r_mid[end])[1]
exact_mid = segment_step_response(t, params_mid.σ, params_mid.z, params_mid.D, params_mid.H, α, kg)
pred_mid - exact_mid

Plots.plot(r_mid, h_midpoint.(r_mid))
####

## Comparison of two approximations
diff(r) = h_mean_sts(r) - h_midpoint(r)
Plots.plot(r_mean, diff.(r_mean))

quadgk(r -> diff(r) * point_step_response(t, r*rb, α, kg), r_mean[1]+10^-14, r_mean[end])[1]
