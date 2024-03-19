
# Solution form the discrete convolution for the load q
function convolve_step(q; Δt, r, α = 10^-6, kg = 3.)
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = point_step_response.(t, r, α, kg)
    return conv(q, response)[1:length(q)]
end

# Step response of the point source
point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)


function convolve_fls_step(q; Δt, r, z, D, H, α = 10^-6, kg = 3.)
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = segment_step_response.(t, r, z, D, H, α, kg)
    return conv(q, response)[1:length(q)]
end

segment_step_response(t, r, z, D, H, α, kg) = quadgk(ζ -> point_step_response(t, sqrt(r^2 + (z-ζ)^2), α, kg), D, D+H)[1]


function convolve_sls_step(q; Δt, r, D1, H1, D2, H2, α = 10^-6, kg = 3.)
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = segment_mean_step_response.(t, r, D1, H1, D2, H2, α, kg)
    return conv(q, response)[1:length(q)]
end

segment_mean_step_response(t, r, D1, H1, D2, H2, α, kg) = hcubature(ζ -> point_step_response(t, sqrt(r^2 + (ζ[1]-ζ[2])^2), α, kg), [D1, D2], [D1+H1, D2+H2])[1] / H2
