include("definitions.jl")
using SpecialFunctions

r(σ, ξ1, ξ2, a1, b1, a2, b2) = sqrt(σ^2 + 1/4 * ( (b1-a1)*ξ1 + (b1+a1) - (b2-a2)*ξ2 - (b2+a2) )^2)
point_response(r, t, α, kg) = erfc(r/sqrt(4α*t)) / (4π*kg*r)
function g_Q_T!(I, bh_disc::BoreholeDiscretization, basis::LagrangeBasis, σ, params::Constants)
    @unpack x, N = basis
    @unpack segments, S, L, D, H = bh_disc
    @unpack α, kg, Δt = params

    sa = @. (segments[1:end-1]+1) * H/2 + D
    sb = @. (segments[2:end]+1) * H/2 + D
    Np = size(I)[1]
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            I[i, j] = quadgk(ξ -> L[v] / 2 * ψ(ξ, n, basis) * point_response(r(σ, ξ, x[m], sa[v], sb[v], sa[u], sb[u]), Δt, α, kg), -1., 1.)[1]
        end
    end
    nothing
end

function g_Tin_Q!(I, params::InternalModelParams, bh_disc::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack segments = bh_disc
    @unpack Rp = params
    Cf = (f1(1, params) + f2(1, params))/(f3(1, params)-f2(1, params))
    Np = length(I)
    for i in 1:Np
        (u, m) = indices(i, N)
        ξ = ξseg(x[m], u, segments)
        I[i] = -Rp/2 * (f1(ξ, params) - f2(ξ, params) +  Cf * (f2(ξ, params) + f3(ξ, params)))
    end
end

function g_T_Q!(I, params::InternalModelParams, bh_disc::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack L, segments = bh_disc
    Np = size(I)[1]
    Cf = 1 / (f3(1., params)-f2(1., params))
    R = params.Rp / 2
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            ξu = ξseg(x[m], u, segments)
            ηv(η) = ξseg(η, v, segments)
            if v < u
                int = quadgk(η -> (Cf * (f2(ξu, params) + f3(ξu, params)) * (f4(-ηv(η), params) + f5(-ηv(η), params)) + f4(ξu - 1 - ηv(η), params) - f5(ξu - 1 - ηv(η), params)) * ψ(η, n, basis), -1., 1.)[1]
                I[i, j] = -R * L[v]/2 * int
            elseif u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                int = quadgk(η -> Cf * (f2(ξu, params) + f3(ξu, params)) * (f4(-ηv(η), params) + f5(-ηv(η), params)) * ψ(η, n, basis) + (x[m]+1)/2 * (f4(ξu - 1 - ηv(ϕ(η)), params) - f5(ξu - 1 - ηv(ϕ(η)), params)) * ψ(ϕ(η), n, basis), -1., 1.)[1]
                I[i, j] = (n == m ? 2R : 0) - R * L[v]/2 * int
            else
                int = quadgk(η -> (Cf * (f2(ξu, params) + f3(ξu, params)) * (f4(-ηv(η), params) + f5(-ηv(η), params))) * ψ(η, n, basis), -1., 1.)[1]
                I[i, j] = -R * L[v]/2 * int
            end
        end
    end
    nothing
end

function g_Q!(I, bh_disc::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack L = bh_disc
    @unpack w = basis
    Np = length(I)
    N = length(w)
    for i in 1:Np
        (v, n) = indices(i, N)
        I[i] = L[v] / 2 * w[n] 
    end
    nothing
end
