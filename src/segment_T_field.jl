using QuadGK

function T_ls_integrand(x, r, z, s::BoreholeSegment) 
    positive_contributions = [z - s.D, z + s.D + s.L]
    negative_contributions = [z + s.D, z - s.D - s.L]

    return exp(- r^2 * x^2) / x * (reduce(+, erf.(positive_contributions * x)) - reduce(+, erf.(negative_contributions * x)))
end

T_ls(r, z, t, s::BoreholeSegment, α = 3*10^-6) = quadgk(x -> T_ls_integrand(x, r, z, s), 1/sqrt(4α * t), Inf)[1]
