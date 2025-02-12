include("definitions.jl")
using SpecialFunctions

r(x, y) = norm(x-y)
#r(σ, ξ1, ξ2, a1, b1, a2, b2) = sqrt(σ^2 + 1/4 * ( (b1-a1)*ξ1 + (b1+a1) - (b2-a2)*ξ2 - (b2+a2) )^2)
point_response(r, t, α, kg) = erfc(r/sqrt(4α*t)) / (4π*kg*r)
function g_Q_T!(I, source::BoreholeDiscretization, target::BoreholeDiscretization, basis::LagrangeBasis, params::Constants, t)
    @unpack x, N = basis
    @unpack segments, S = bh_disc
    @unpack α, kg, rb = params

    sr = source == target ? [rb, 0., 0.] :  [0., 0., 0.]

    Np = size(I)[1]
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N) # target
            (v, n) = indices(j, N) # source
            #I[i, j] = quadgk(ξ -> normJ(ξ, v, source) * ψ(ξ, n, basis) * point_response(r(σ, ξ, x[m], sa[v], sb[v], sa[u], sb[u]), t, α, kg), -1., 1.)[1]
            I[i, j] = quadgk(ξ -> normJ(ξ, v, source) * ψ(ξ, n, basis) * point_response(r(path(ξ, v, source) + sr, path(x[m], u, target)), t, α, kg), -1., 1.)[1]
        end
    end
    nothing
end

function g_Tin_Q!(I, params::InternalModelParams, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack segments = bh
    @unpack Rp = params
    Cf = (f1(1, params, bh) + f2(1, params, bh))/(f3(1, params, bh)-f2(1, params, bh))
    Np = length(I)
    for i in 1:Np
        (u, m) = indices(i, N)
        ξ = ξseg(x[m], u, segments)
        I[i] = -Rp/2 * (f1(ξ, params, bh) - f2(ξ, params, bh) +  Cf * (f2(ξ, params, bh) + f3(ξ, params, bh)))
    end
end

function g_T_Q!(I, params::InternalModelParams, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack segments = bh
    Np = size(I)[1]
    Cf = 1 / (f3(1., params, bh)-f2(1., params, bh))
    R = params.Rp / 2
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            ξu = ξseg(x[m], u, segments)
            ηv(η) = ξseg(η, v, segments)
            if v < u
                int = quadgk(η -> (Cf * (f2(ξu, params, bh) + f3(ξu, params, bh)) * (f4(-ηv(η), params, bh) + f5(-ηv(η), params, bh)) + f4(ξu - 1 - ηv(η), params, bh) - f5(ξu - 1 - ηv(η), params, bh)) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
                I[i, j] = -R * int
            elseif u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                int = quadgk(η -> Cf * (f2(ξu, params, bh) + f3(ξu, params, bh)) * (f4(-ηv(η), params, bh) + f5(-ηv(η), params, bh)) * ψ(η, n, basis) * normJ(η, v, bh) + (x[m]+1)/2 * (f4(ξu - 1 - ηv(ϕ(η)), params, bh) - f5(ξu - 1 - ηv(ϕ(η)), params, bh)) * ψ(ϕ(η), n, basis) * normJ(ϕ(η), v, bh), -1., 1.)[1]
                I[i, j] = (n == m ? 2R : 0) - R  * int
            else
                int = quadgk(η -> (Cf * (f2(ξu, params, bh) + f3(ξu, params, bh)) * (f4(-ηv(η), params, bh) + f5(-ηv(η), params, bh))) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
                I[i, j] = -R * int
            end
        end
    end
    nothing
end

function g_Q!(I, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, w = basis
    Np = length(I)
    N = length(w)
    for i in 1:Np
        (v, n) = indices(i, N)
        I[i] = normJ(x[n], v, bh) *  w[n] 
    end
    nothing
end
