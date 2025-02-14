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

function g_Tin_Q!(I, internal::InternalModelParams, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack segments = bh
    #@unpack Rp = params
    @unpack Rb, R11Δ, R22Δ = internal
    Cf = (f1(1, internal, bh) + f2(1, internal, bh))/(f3(1, internal, bh)-f2(1, internal, bh))
    Np = length(I)
    for i in 1:Np
        (u, m) = indices(i, N)
        ξ = ξseg(x[m], u, segments)
        #I[i] = -Rp/2 * (f1(ξ, internal, bh) - f2(ξ, internal, bh) +  Cf * (f2(ξ, internal, bh) + f3(ξ, internal, bh)))
        I[i] = - ( (f1(ξ, internal, bh) + Cf * f2(ξ, internal, bh))/R11Δ + (-f2(ξ, internal, bh) + Cf * f3(ξ, internal, bh))/R22Δ )
    end
end

function g_T_Q!(I, internal::InternalModelParams, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, N = basis
    @unpack segments = bh
    @unpack Rb, R11Δ, R22Δ = internal
    Np = size(I)[1]
    C = 1 / (f3(1., internal, bh) - f2(1., internal, bh))
    #R = internal.Rp / 2
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            ξu = ξseg(x[m], u, segments)
            ηv(η) = ξseg(η, v, segments)
            #=if v < u
                int = quadgk(η -> (Cf * (f2(ξu, internal, bh) + f3(ξu, internal, bh)) * (f4(-ηv(η), internal, bh) + f5(-ηv(η), internal, bh)) + f4(ξu - 1 - ηv(η), internal, bh) - f5(ξu - 1 - ηv(η), internal, bh)) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
                I[i, j] = -R * int
            elseif u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                int = quadgk(η -> Cf * (f2(ξu, internal, bh) + f3(ξu, internal, bh)) * (f4(-ηv(η), internal, bh) + f5(-ηv(η), internal, bh)) * ψ(η, n, basis) * normJ(η, v, bh) + (x[m]+1)/2 * (f4(ξu - 1 - ηv(ϕ(η)), internal, bh) - f5(ξu - 1 - ηv(ϕ(η)), internal, bh)) * ψ(ϕ(η), n, basis) * normJ(ϕ(η), v, bh), -1., 1.)[1]
                I[i, j] = (n == m ? 2R : 0) - R  * int
            else
                int = quadgk(η -> (Cf * (f2(ξu, internal, bh) + f3(ξu, internal, bh)) * (f4(-ηv(η), internal, bh) + f5(-ηv(η), internal, bh))) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
                I[i, j] = -R * int
            end=#

            F(η) = - C * ( f2(ξu, internal, bh) / R11Δ + f3(ξu, internal, bh) / R22Δ ) * ( f4(-ηv(η), internal, bh) + f5(-ηv(η), internal, bh) ) * ψ(η, n, basis) * normJ(η, v, bh)
            if u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                F_self_s(η) = - (x[m]+1)/2 * (f4(ξu - 1 - ηv(ϕ(η)), internal, bh) / R11Δ - f5(ξu - 1 - ηv(ϕ(η)), internal, bh) / R22Δ) * ψ(ϕ(η), n, basis) * normJ(ϕ(η), v, bh)
                I[i, j] = (n == m ? 1/Rb : 0) + quadgk(η -> F(η) + F_self_s(η), -1., 1.)[1]
            elseif v < u
                F_sts(η) = - (f4(ξu - 1 - ηv(η), internal, bh) / R11Δ - f5(ξu - 1 - ηv(η), internal, bh) / R22Δ) * ψ(η, n, basis) * normJ(η, v, bh)
                I[i, j] = quadgk(η -> F(η) + F_sts(η), -1., 1.)[1]
            else 
                I[i, j] = quadgk(F, -1., 1.)[1]
            end
        end
    end
    nothing
end

function g_Q!(I, bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack x, w = basis
    Np = length(I)
    N = length(w)
    L = s(1., bh)
    for i in 1:Np
        (v, n) = indices(i, N)
        I[i] = normJ(x[n], v, bh) *  w[n] / L 
    end
    nothing
end


function g_out!(I, basis::LagrangeBasis, bh::BoreholeDiscretization, params::InternalModelParams) 
    @unpack x, N = basis
    @unpack segments = bh
    Np = size(I)[1]

    Cf = 1 / (f3(1., params, bh)-f2(1., params, bh))
    for i in 1:Np
        (v, n) = indices(i, N)
        I[i] = Cf * quadgk(η -> (f4(-ξseg(η, v, segments), params, bh) + f5(-ξseg(η, v, segments), params, bh)) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1] 
    end
end

function g_T!(I, f, basis::LagrangeBasis, bh::BoreholeDiscretization) 
    @unpack x, N = basis
    @unpack segments, S = bh
    Np = N*S
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            ξu = ξseg(x[m], u, segments)
            ηv(η) = ξseg(η, v, segments)
            if v < u
                I[i, j] = quadgk(η -> f(ξu - 1 - ηv(η)) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
            elseif u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                I[i, j] = (x[m]+1)/2 * quadgk(η -> f(ξu - 1 - ηv(ϕ(η))) * ψ(ϕ(η), n, basis) * normJ(ϕ(η), v, bh), -1., 1.)[1]
            else 
                I[i, j] = 0.
            end
        end
    end
end
g_T1!(I, basis::LagrangeBasis, bh::BoreholeDiscretization, params::InternalModelParams) = g_T!(I, x -> f4(x, params, bh), basis, bh) 
g_T2!(I, basis::LagrangeBasis, bh::BoreholeDiscretization, params::InternalModelParams) = g_T!(I, x -> f5(x, params, bh), basis, bh) 

function fluid_profiles(Tin, Tb, bh::BoreholeDiscretization, internal::InternalModelParams, basis::LagrangeBasis)
    Np = bh.S * basis.N
    G_Tin_Tout = (f1(1., internal, bh)+f2(1., internal, bh)) / (f3(1., internal, bh)-f2(1., internal, bh))
    Gout = zeros(Np)
    g_out!(Gout, basis, bh, internal)
    Tout = G_Tin_Tout * Tin + dot(Gout, Tb)

    GT1 = zeros(Np, Np)
    GT2 = zeros(Np, Np)
    g_T1!(GT1, basis, bh, internal) 
    g_T2!(GT2, basis, bh, internal) 

    ξ_disc = ξ_discretization(bh, basis)

    f1_ξ = f1.(ξ_disc, Ref(internal), Ref(bh))
    f2_ξ = f2.(ξ_disc, Ref(internal), Ref(bh))
    f3_ξ = f3.(ξ_disc, Ref(internal), Ref(bh))

    T1 =  f1_ξ .* Tin + f2_ξ .* Tout + GT1*Tb
    T2 = -f2_ξ .* Tin + f3_ξ .* Tout - GT2*Tb

    (Tout, T1, T2)
end
