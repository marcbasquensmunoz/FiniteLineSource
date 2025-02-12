include("definitions.jl")
include("coefficients.jl")

N = 5
H = 100.
D = 0.
segments = [-1., -0.6, 0.6, 1.]#[-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization([-1., 1.], [[0., 0., D], [0., 0., D+H]], segments) 

basis = LagrangeBasis(N)
Np = basis.N * bh_disc.S

mf = 0.05
cpf = 4000.
R11 = 0.6
R22 = 0.3
R12 = -5.

params = InternalModelParams(mf=mf, cpf=cpf, R11=R11, R22=R22, R12=R12)

ξ_disc = reduce(vcat, [ξseg.(basis.x, Ref(u), Ref(segments)) for u in 1:bh_disc.S])

Rp = (R11+R22)/(R11*R22)

β1 = 1/(mf*cpf*R11)
β2 = 1/(mf*cpf*R22)
β12 = 1/(mf*cpf*R12)

β = (β2-β1)/2
γ = sqrt((β2+β1)^2 / 4 + (β2+β1)*β12)
δ = 1/γ * (β12 + (β2+β1)/2)
s_ξ = @. H/2 * (ξ_disc+1) + D
f1_ξ = @. exp(β*s_ξ) * (cosh(γ*s_ξ) - δ*sinh(γ*s_ξ))
f2_ξ = @. exp(β*s_ξ) * β12/γ * sinh(γ*s_ξ)
f3_ξ = @. exp(β*s_ξ) * (cosh(γ*s_ξ) + δ*sinh(γ*s_ξ))
f4_ξ = @. exp(β*s_ξ) * (β1*cosh(γ*s_ξ) - (δ*β1 + β2*β12/γ)*sinh(γ*s_ξ))
f5_ξ = @. exp(β*s_ξ) * (β2*cosh(γ*s_ξ) + (δ*β2 + β1*β12/γ)*sinh(γ*s_ξ))

G_Tin_Tout = (f1(1., params, bh_disc)+f2(1., params, bh_disc)) / (f3(1., params, bh_disc)-f2(1., params, bh_disc))

##################
# Fluid profile
##################

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
    @unpack segments = bh
    Np = size(I)[1]
    for j in 1:Np
        for i in 1:Np
            (u, m) = indices(i, N)
            (v, n) = indices(j, N)
            ξu = ξseg(x[m], u, segments)
            if v < u
                I[i, j] = quadgk(η -> f(ξu - 1 - ξseg(η, v, segments)) * ψ(η, n, basis) * normJ(η, v, bh), -1., 1.)[1]
            elseif u == v
                ϕ(η) = (x[m]+1)/2 * η + (x[m]-1)/2
                I[i, j] = quadgk(η -> (x[m]+1)/2 * f(ξu - 1 - ξseg(ϕ(η), v, segments)) * ψ(ϕ(η), n, basis) * normJ(ϕ(η), v, bh), -1., 1.)[1]
            else 
                I[i, j] = 0.
            end
        end
    end
end
g_T1!(I, basis::LagrangeBasis, bh::BoreholeDiscretization, params::InternalModelParams) = g_T!(I, x -> f4(x, params, bh), basis, bh) 
g_T2!(I, basis::LagrangeBasis, bh::BoreholeDiscretization, params::InternalModelParams) = g_T!(I, x -> f5(x, params, bh), basis, bh) 

Gout = zeros(Np)
GT1 = zeros(Np, Np)
GT2 = zeros(Np, Np)

g_out!(Gout, basis, bh_disc, params)
g_T1!(GT1, basis, bh_disc, params) 
g_T2!(GT2, basis, bh_disc, params) 


Tb = vcat(8*ones(N*1), 10*ones(2*N))
Tin = 10.
Tout = G_Tin_Tout*Tin + dot(Gout, Tb) 
T1 =  f1_ξ .* Tin + f2_ξ .* Tout + GT1*Tb
T2 = -f2_ξ .* Tin + f3_ξ .* Tout - GT2*Tb
q = (Tb - 1/2*(T1+T2)) * Rp

z = @. H/2 * (ξ_disc+1) + D
scatter(T1, z, label="T1", yflip = true, xlimits=(0,15.))
scatter!(T2, z, label="T2")
scatter!(q, z, label="q")
plot!(Tb, z, label="Tb")


###################
# Coefficient check
###################

GTinQ = zeros(Np)
GTQ = zeros(Np, Np)
g_Tin_Q!(GTinQ, params, bh_disc, basis)
g_T_Q!(GTQ, params, bh_disc, basis)

check_Tin = - (f1_ξ - f2_ξ + (f2_ξ + f3_ξ)*G_Tin_Tout ) * Rp / 2
check_Tin - GTinQ

check_G = -(GT1-GT2 + (f2_ξ + f3_ξ) * Gout' - 2*Diagonal(ones(Np))) * Rp / 2
check_G - GTQ


########################
# Heat extraction from T
########################
Tin = 10.
Tb = ones(Np)
q = GTinQ .* Tin + GTQ * Tb

plot(ξ_disc, q)