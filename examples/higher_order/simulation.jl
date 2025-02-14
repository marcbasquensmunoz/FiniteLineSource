include("definitions.jl")
include("coefficients.jl")
using Plots

α = 1e-6
kg = 3.
rb = 0.1
Δt = 30*24*3600.
Nt = 100
T0 = 10.

constants = Constants(α=α, kg=kg, rb=rb, Δt=Δt)

D = 0.
H = 100.
segments = [-1., -0.6, 0.6, 1.]#[-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization([-1., 1.], [[0., 0., D], [0., 0., D+H]], segments) 

mf = 0.05
cpf = 4000.
R11Δ = 0.6
R22Δ = 0.3
R12Δ = -5.
internal = InternalModelParams(mf=mf, cpf=cpf, R11Δ=R11Δ, R22Δ=R22Δ, R12Δ=R12Δ)

N = 5
basis = LagrangeBasis(N)

Np = N*bh_disc.S

ξ_disc = ξ_discretization(bh_disc, basis)

Qtot = -1.

# Allocate
Ht = zeros(Np, Np, Nt)
M = zeros(2Np+1, 2Np+1)
b = zeros(2Np+1)
X = zeros(2Np+1, Nt)

rTb = 1:Np
rq = Np+1:2Np
rTin = 2Np+1

# Precompute ground response
for i in 1:Nt
    @views g_Q_T!(Ht[:, :, i], bh_disc, bh_disc, basis, constants, i*Δt)
end 

# Fill matrix
for i in 1:Np
    M[i, i] = 1. 
    M[Np+i, Np+i] = -1. 
end
@views @. M[rTb, rq] = Ht[:, :, 1]
@views g_T_Q!(M[rq, rTb], internal, bh_disc, basis)
@views g_Tin_Q!(M[rq, rTin], internal, bh_disc, basis)
@views g_Q!(M[rTin, rq], bh_disc, basis)
b[rTb] .= T0
b[rTin] = Qtot

for i in 1:Nt
    # Update 
    b[rTb] .= T0
    for j in 1:i-1
        @views b[rTb] .-= (Ht[:, :, i-j+1] - Ht[:, :, i-j]) * X[rq, j]
    end
    
    # Solve each time step
    X[:, i] .= M \ b
end

Tb = X[rTb, end]
q = X[rq, end]
Tin = X[rTin, end]
Tout, T1, T2 = fluid_profiles(Tin, Tb, bh_disc, internal, basis)

Tbz(ξ) = evaluate(ξ, Tb, bh_disc, basis)
qz(ξ) = evaluate(ξ, q, bh_disc, basis)

# Fluid profiles
scatter(ξ_disc, T1)
scatter!(ξ_disc, T2)
plot!(ξ_disc, Tb)

# Heat profile
ξξ = -1:0.01:1
plot(ξξ, qz.(ξξ))

# Energy balance
Q_int = integrate_bh(q, bh_disc, basis)
Q_balance = (Tout-Tin)*internal.cpf*internal.mf / H
