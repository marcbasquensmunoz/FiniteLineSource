
# Solution form the discrete convolution for the load q
function convolve_step(q; Δt, r, α = 10^-6, kg = 3.)
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = point_step_response.(t, r, α, kg)
    return conv(q, response)[1:length(q)]
end

# Step response of the point source
point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)


function convolve_step(q, model::Setup; params::Constants)
    @unpack Δt = params
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = step_response.(t, Ref(model), Ref(params))
    return conv(q, response)[1:length(q)]
end

# Point to point
function step_response(t, model::PointToPoint, params::Constants)
    @unpack r = model
    @unpack α, kg = params
    point_step_response(t, r, α, kg)
end

# Segment to point
function step_response(t, model::SegmentToPoint, params::Constants)
    @unpack r, D, H, z = model
    @unpack α, kg = params
    quadgk(ζ -> point_step_response(t, sqrt(r^2 + (z-ζ)^2), α, kg), D, D+H)[1]
end

# Mean segment to segment
function step_response(t, model::SegmentToSegment, params::Constants)
    @unpack r, D1, H1, D2, H2 = model
    @unpack α, kg = params
    params = MeanSegToSegEvParams(model)
    h_mean_sts, r_min, r_max = mean_sts_evaluation(params)
    quadgk(r -> h_mean_sts(r) * point_step_response(t, r, α, kg), r_min, r_max)[1]
end

# Moving point source
function step_response(t, model::MovingPointToPoint, params::Constants)
    @unpack r, x, v = model
    @unpack α, kg = params
    exp(v * (x-r)/ (2α)) * (erfc( (r-t*v) / sqrt(4t*α)) + exp(v*r/α) * erfc((r+t*v) / sqrt(4t*α)) ) / (8π*r*kg)
end