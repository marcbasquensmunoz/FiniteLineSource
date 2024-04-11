
function compute_self_response(n; t, rb, H, α, kg)
    σ = rb
    params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=σ)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)

    principal(r) = erfc(σ/sqrt(4*t*α)) / (2π*kg*sqrt(2*σ*(r-σ)))
    f(r) = h_mean_sts(r) * point_step_response(t, r, α, kg) - principal(r)
    principal_integral = sqrt(r_max - r_min) * erfc(r_min/sqrt(4α*t)) / (kg*π*sqrt(2*r_min))

    integrate(f, r_min, r_max, n) + principal_integral
end
