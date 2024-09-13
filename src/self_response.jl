
function compute_self_response(n; t, rb, H, α, kg)
    σ = rb

    x, w = gausslegendre(n+1)

    params = MeanSegToSegEvParams(D1=0., H1=H, D2=0., H2=H, σ=σ)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)

    principal(r) = erfc(σ/sqrt(4*t*α)) / (2π*kg*sqrt(2*σ*(r-σ)))
    f(r) = h_mean_sts(r) * point_step_response(t, r, α, kg) - principal(r)
    principal_integral = sqrt(r_max - r_min) * erfc(r_min/sqrt(4α*t)) / (kg*π*sqrt(2*r_min))
    b = 5*rb

    integrate(f, r_min, b, x, w) + integrate(f, b, r_max, x, w) + principal_integral
    #integrate(f, r_min, b, n) + integrate(f, b, r_max, n) + principal_integral
end

function integrate(f, a, b, x, w)
    m = (b-a)/2
    c = (b+a)/2
    m * dot(w, @. f(m*x+c))
end

function compute_sts_response(s:: SegmentToSegment, params::Constants)
    params = MeanSegToSegEvParams(s)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)
    f(r) = h_mean_sts(r) * point_step_response(t, r, α, kg)
    x, w = adaptive_gk(f, r_min, r_max)
    dot(f.(x), w)
end
