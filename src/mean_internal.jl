
@with_kw struct InternalSegToSegEvParams{T <: Number} @deftype T
    D1
    H1
    D2
    H2
    A::Matrix
    σ

    rmin = σ
    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     )
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     )
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin
end
transpose(p::InternalSegToSegEvParams) = InternalSegToSegEvParams(D1=p.D2, H1=p.H2, D2=p.D1, H2=p.H1, σ=p.σ, A=p.A)

function Lnew(r, params::InternalSegToSegEvParams)
    @unpack D1, H1, D2, H2, σ, r1, r2, r3, r4 = params

    if r < r1
        0.
    elseif r < r2
        (D1 - D2 + H1) + sqrt(r^2-σ^2)
    elseif r < r3
        min(H1, H2)
    elseif r < r4
        (D2 - D1 + H2) - sqrt(r^2-σ^2)
    else
        0.
    end
end

function mean_internal_evaluation(params::InternalSegToSegEvParams)
    @unpack A, σ = params
    paramsT = transpose(params)
    r_min = Base.max(params.r1, paramsT.r1)
    r_max = Base.max(params.r4, paramsT.r4)    
    Lt(r) = Lnew(r, params) + Lnew(r, paramsT)
    e = sum([1 -1] * exp(A*params.H2) * [1, 1])
    h_new(r) = sum([1 -1] * r/sqrt(r^2-σ^2) * exp(A*sqrt(r^2-σ^2)) * (I - exp(-A*Lt(r))) * [1, 1]) / e
    h_new, r_min, r_max
end