
@with_kw struct LevelSetParams
    D1
    H1
    D2
    H2
    σ
    rb

    rmin = σ/rb
    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     )/rb 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )/rb 
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     )/rb 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)/rb 

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin
end

transpose(p::LevelSetParams) = LevelSetParams(D1=p.D2, H1=p.H2, D2=p.D1, H2=p.H1, σ=p.σ, rb=p.rb)

grad(r, rb, σ) = sqrt(2)/(rb^2*r)*sqrt((rb*r)^2-σ^2)

function L(r, params::LevelSetParams)
    @unpack D1, H1, D2, H2, σ, rb, r1, r2, r3, r4 = params

    if r < r1
        0
    elseif r < r2
        sqrt(2)*(D1 - D2 + H1 + sqrt((r*rb)^2-σ^2))
    elseif r < r3
        sqrt(2)*min(H1, H2)
    elseif r < r4
        sqrt(2)*(D2 - D1 + H2 - sqrt((r*rb)^2-σ^2))
    else
        0
    end
end

function mean_sts_evaluation(params::LevelSetParams)
    paramsT = transpose(params)
    f(r) = L(r, params) + L(r, paramsT)
    h_mean_sts(r) = f(r) / grad(r, params.rb, params.σ) / params.H2
    r = max(params.r1, paramsT.r1):0.01:max(params.r4, paramsT.r4)
    return h_mean_sts, r    
end
