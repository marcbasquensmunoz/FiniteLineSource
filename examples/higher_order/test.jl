include("definitions.jl")
include("coefficients.jl")
using Plots
using FastChebInterp
using Roots

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
bh_disc = BoreholeDiscretization([-1., 0., 1.], [[0., 0., D], [5., 0., D+H/2], [0., 0., D+H]], segments) 

xi = -1:0.01:1
pathxi = path.(xi, Ref(bh_disc))
X1 = [xx[1] for xx in pathxi]
X3 = [xx[3] for xx in pathxi]

plot(X1, X3)

L(ξ) = quadgk(η -> s(η, bh_disc), -1, ξ)[1]

plot(xi, L.(xi))

t(s0) = find_zero(ξ -> s(ξ, bh_disc) - s0, 0.)

x = chebpoints(200, 0, L(1))
c = chebinterp(t.(x), 0, L(1))

γ(s) = path(c(s), bh_disc)
dγ(s) = norm(dpath(c(s), bh_disc) * chebgradient(c, s)[2])
