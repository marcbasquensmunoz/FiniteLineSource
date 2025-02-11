include("definitions.jl")
include("coefficients.jl")

α = 1e-6
kg = 3.
rb = 0.1

D = 0.
H = 100.
σ = 0.1

t = 30*24*3600.

N = 10
basis = LagrangeBasis(N)
segments = [-1, -0.95, -0.9, -0.6, 0.6, 0.9, 0.95, 1.]
bh_disc = BoreholeDiscretization(D, H, segments) 
constants = Constants(α=α, kg=kg, Δt=t, rb=rb)

Np = bh_disc.S * N

q = vcat(2*ones(N*3), ones(N*4))
GQT = zeros(Np, Np)
g_Q_T!(GQT, bh_disc, basis, σ, constants)
Tk = GQT * q

response(z1, z2) = erfc(sqrt(σ^2 + (z1-z2)^2)/sqrt(4α*t)) / (4π*kg*sqrt(σ^2 + (z1-z2)^2))
Ttest(ξ) = quadgk(η -> H/2 * evaluate(η, q, bh_disc, basis) * response(H/2*(ξ+1)+D, H/2*(η+1)+D) , -1., 1.)[1]

ξξ = -1:0.01:1
T1 = evaluate.(ξξ, Ref(Tk), Ref(bh_disc), Ref(basis)) 
T2 = Ttest.(ξξ)

plot(log10.(abs.(T1-T2)), ylim=(-10,0))