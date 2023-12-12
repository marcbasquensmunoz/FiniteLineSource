using SpecialFunctions
using QuadGK

ierf(x) = x * erf(x) - (1 - exp(-x^2)) / sqrt(π)

"Integrand of the thermal response for the borehole wall of s₂ caused by heat extraction at s₁"
function seg2seg_integrand(x, s₁::BoreholeSegment, s₂::BoreholeSegment) 
    r² = s₁ == s₂ ? s₁.r_b^2 : (s₁.x - s₂.x)^2 + (s₁.y - s₂.y)^2 
    positive_contributions = [s₂.D - s₁.D + s₂.L, s₂.D - s₁.D - s₁.L, s₁.D + s₂.D + s₁.L,  s₁.D + s₂.D + s₂.L]
    negative_contributions = [s₂.D + s₁.D, s₂.D - s₁.D, s₂.D - s₁.D + s₂.L - s₁.L, s₂.D + s₁.D + s₂.L + s₁.L]
    return exp(-r² * x^2) / x^2 * (reduce(+, ierf.(positive_contributions * x)) - reduce(+, ierf.(negative_contributions * x)))
end

function segment_to_segment_step_response(t; source::BoreholeSegment, receptor::BoreholeSegment, α = 3*10^-6)
    t == 0 ? 0.0 : 1 / (2 * receptor.L) * quadgk(x -> seg2seg_integrand(x, source, receptor), 1/sqrt(4α * t), Inf)[1]
end
