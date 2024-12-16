
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
    @unpack σ, D, H, z = model
    @unpack α, kg = params

    I_stp(s) = erf(s * (z-D)) - erf(s * (z-D-H))
    1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * I_stp(s), 1/sqrt(4*α*t), Inf, atol = eps())[1]

#    quadgk(ζ -> point_step_response(t, sqrt(σ^2 + (z-ζ)^2), α, kg), D, D+H, atol=eps())[1]
end

# Mean segment to segment
function step_response(t, model::SegmentToSegment, params::Constants)
    @unpack σ, D1, H1, D2, H2 = model
    @unpack α, kg, rb = params
    params = MeanSegToSegEvParams(model)
    r_min, r_max = h_mean_lims(params)
    quadgk(r -> h_mean_sts(r, params) * point_step_response(t, r, α, kg), r_min, r_max, atol=eps())[1]
end

function step_response(t, model::SegmentToSegmentOld, params::Constants)
    @unpack σ, D1, H1, D2, H2 = model
    @unpack α, kg = params

    I(s) = ierf((D2 - D1 + H2)*s) + ierf((D2 - D1 - H1)*s) - ierf((D2 - D1)*s) - ierf((D2 - D1 + H2 - H1)*s) 
    return 1/(4π*kg) * quadgk(s -> exp(-σ^2 * s^2) / s^2 * I(s), 1/sqrt(4α*t), Inf, atol=eps())[1]

    #hcubature(z -> point_step_response(t, sqrt(σ^2 + (z[1]-z[2])^2), α, kg), [D1, D2], [D1+H1, D2+H2], abstol=1e-10)[1]
    atol = eps()
    quadgk(z1 -> quadgk(z2 -> point_step_response(t, sqrt(σ^2 + (z1-z2)^2), α, kg), D2, D2+H2, atol=atol)[1], D1, D1+H1, atol=atol)[1]
end 

# Moving point source
function step_response(t, model::MovingPointToPoint, params::Constants)
    @unpack r, x, v = model
    @unpack α, kg = params
    if v*r/α > 700 && (r+t*v) / sqrt(4t*α) > 6 # Avoid NaN
        return exp(-v * (r-x)/ (2α)) * (erfc( (r-t*v) / sqrt(4t*α))) / (8π*r*kg)
    end
    exp(-v * (r-x)/ (2α)) * (erfc( (r-t*v) / sqrt(4t*α)) + erfc((r+t*v) / sqrt(4t*α)) * exp(v*r/α) ) / (8π*r*kg)
end

function step_response(t, model::MovingSegmentToPoint, params::Constants)
    @unpack σ, x, z, v, D, H = model
    @unpack α, kg = params
    function f(ζ) 
        val = exp(v * (x-sqrt(σ^2 + (z-ζ)^2))/ (2α)) * (erfc( (sqrt(σ^2 + (z-ζ)^2)-t*v) / sqrt(4t*α)) + exp(v*sqrt(σ^2 + (z-ζ)^2)/α) * erfc((sqrt(σ^2 + (z-ζ)^2)+t*v) / sqrt(4t*α)) ) / (8π*sqrt(σ^2 + (z-ζ)^2)*kg)
        if isnan(val)
            return 0.
        end
        val
    end
    quadgk(ζ -> f(ζ), D, D+H)[1]
end

ierf(x) = x*erf(x) - (1 - exp(-x^2)) / sqrt(π)