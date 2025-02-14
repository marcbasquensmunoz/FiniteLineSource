include("definitions.jl")
include("coefficients.jl")
using Plots

N = 5
H = 100.
D = 0.
segments = [-1., -0.6, 0.6, 1.]#[-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization([-1., 1.], [[0., 0., D], [0., 0., D+H]], segments) 

basis = LagrangeBasis(N)
Np = basis.N * bh_disc.S

mf = 0.05
cpf = 4000.
R11Δ = 0.6
R22Δ = 0.3
R12Δ = -5.

internal = InternalModelParams(mf=mf, cpf=cpf, R11Δ=R11Δ, R22Δ=R22Δ, R12Δ=R12Δ)

ξ_disc = ξ_discretization(bh_disc, basis)

Tb = vcat(8*ones(N*1), 10*ones(2*N))
Tin = 10.
Tout, T1, T2 = fluid_profiles(Tin, Tb, bh_disc, internal, basis)
q = (Tb - 1/2*(T1+T2)) * internal.Rp

z = @. H/2 * (ξ_disc+1) + D
scatter(T1, z, label="T1", yflip = true, xlimits=(0,15.))
scatter!(T2, z, label="T2")
plot!(Tb, z, label="Tb")

scatter(q, z, label="q")

Q_balance = (Tout-Tin)*cpf*mf 
Q_int = integrate_bh(q, bh_disc, basis)


###################
# Coefficient check
###################

f1_ξ = f1.(ξ_disc, Ref(internal), Ref(bh_disc))
f2_ξ = f2.(ξ_disc, Ref(internal), Ref(bh_disc))
f3_ξ = f3.(ξ_disc, Ref(internal), Ref(bh_disc))

GTinQ = zeros(Np)
GTQ = zeros(Np, Np)
GT1 = zeros(Np, Np)
GT2 = zeros(Np, Np)
Gout = zeros(Np)

G_Tin_Tout = (f1(1., internal, bh_disc)+f2(1., internal, bh_disc)) / (f3(1., internal, bh_disc)-f2(1., internal, bh_disc))
g_Tin_Q!(GTinQ, internal, bh_disc, basis)
g_T_Q!(GTQ, internal, bh_disc, basis)
g_out!(Gout, basis, bh_disc, internal)
g_T1!(GT1, basis, bh_disc, internal) 
g_T2!(GT2, basis, bh_disc, internal) 


check_Tin = - (f1_ξ - f2_ξ + (f2_ξ + f3_ξ)*G_Tin_Tout ) * internal.Rp / 2
check_G = -(GT1-GT2 + (f2_ξ + f3_ξ) * Gout' - 2*Diagonal(ones(Np))) * internal.Rp / 2

@show sum(abs.(check_Tin - GTinQ))
@show sum(abs.(check_G - GTQ))


########################
# Heat extraction from T
########################
Tin = 10.
Tb = ones(Np)
q = GTinQ .* Tin + GTQ * Tb

plot(ξ_disc, q)