include("definitions.jl")
include("coefficients.jl")
using Plots

α = 1e-6
kg = 3.
rb = 0.1
Δt = 30*24*3600.
T0 = 10.

constants = Constants(α=α, kg=kg, rb=rb, Δt=Δt)

D = 0.
H = 100.
σ = 0.1

mf = 0.05
cpf = 4000.
R11 = 0.6
R22 = 0.3
R12 = -5.
internal = InternalModelParams(mf=mf, cpf=cpf, R11=R11, R22=R22, R12=R12)

N = 10
basis = LagrangeBasis(N)

segments = [-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization(D, H, segments) 
Np = N*bh_disc.S

ξ_disc = reduce(vcat, [ξseg.(x, Ref(u), Ref(segments)) for u in 1:bh_disc.S])

Qtot = H


M = zeros(2Np+1, 2Np+1)
b = zeros(2Np+1)

rTb = 1:Np
rq = Np+1:2Np
rTin = 2Np+1

for i in 1:Np
    M[i, i] = -1. 
    M[Np+i, Np+i] = -1. 
end
@views g_Q_T!(M[rTb, rq], bh_disc, basis, σ, constants)
@views g_T_Q!(M[rq, rTb], internal, bh_disc, basis)
@views g_Tin_Q!(M[rq, rTin], internal, bh_disc, basis)
@views g_Q!(M[rTin, rq], bh_disc, basis)
b[rTb] .= -T0
b[rTin] = Qtot

X = M \ b
Tb = X[1:Np]
q = X[Np+1:2Np]
Tin = X[end]

scatter(ξ_disc, Tb, label = "Tb")
scatter!(ξ_disc, q, label = "q")
