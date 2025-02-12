include("definitions.jl")
include("coefficients.jl")
using Plots

α = 1e-6
kg = 3.
rb = 0.1
Δt = 30*24*3600.
Nt = 12
T0 = 10.

constants = Constants(α=α, kg=kg, rb=rb, Δt=Δt)

D = 0.
H = 100.
segments = [-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization([-1., 1.], [[0., 0., D], [0., 0., D+H]], segments) 

mf = 0.05
cpf = 4000.
R11 = 0.6
R22 = 0.3
R12 = -5.
internal = InternalModelParams(mf=mf, cpf=cpf, R11=R11, R22=R22, R12=R12)

N = 5
basis = LagrangeBasis(N)

Np = N*bh_disc.S

ξ_disc = reduce(vcat, [ξseg.(basis.x, Ref(u), Ref(segments)) for u in 1:bh_disc.S])

Qtot = H

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

Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=1 f5(-0., internal, bh_disc)
PProf.Allocs.pprof(from_c=false)

# Fill matrix
for i in 1:Np
    M[i, i] = 1. 
    M[Np+i, Np+i] = -1. 
end
@views @. M[rTb, rq] = -Ht[:, :, 1]
@views g_T_Q!(M[rq, rTb], internal, bh_disc, basis)
@views g_Tin_Q!(M[rq, rTin], internal, bh_disc, basis)
@views g_Q!(M[rTin, rq], bh_disc, basis)
b[rTb] .= T0
b[rTin] = Qtot

for i in 1:Nt
    # Solve each time step
    X[:, i] .= M \ b
    # Update 
    b[rTb] .= T0
    for j in 1:i-1
        @views b[rTb] .+= (Ht[:, :, i-j+1] - Ht[:, :, i-j]) * X[Np+1:2Np, j]
    end
end


Tb = X[1:Np, end]
q = X[Np+1:2Np, end]
Tin = X[end, end]
Tbz(z) = evaluate(z, Tb, bh_disc, basis)
qz(z) = evaluate(2/H*(z-D)-1, q, bh_disc, basis)

zz = -1:0.01:1
plot(zz, Tbz.(zz))

scatter(ξ_disc, Tb, label = "Tb")
scatter!(ξ_disc, q, label = "q")
