using FiniteLineSource
using FiniteLineSource: PointSource
using BenchmarkTools
using SparseArrays
using FMM3D

ϵ = 1e-6
Δt = 3600.
Nt = 8760

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

B = 5.

positions = [PointSource(B*i, B*j, B*k, rb) for i in 0:1 for j in 0:0 for k in 0:0]
q = ones(length(positions)) * [5 * sin(i/24) for i in 1:Nt]' #ones(length(positions), Nt)

#####################################
# Error analysis
#####################################
setup = PointToPoint(r = B)

# Block method
sources = targets = hcat([[p.x,p.y,p.z] for p in  positions]...)
res = zeros(size(targets)[2])

block = prepare_containers(setup, positions, ϵ, Nt, constants, compute_first_block=false);
Ib = zeros(length(positions), Nt)
evolve!(Ib, q, block)


bound(N1, N2, m) = 1/2 * (log((N2+N1*N2)/(N1+N1*N2)) - 1/m * log((N2*(N1+m))/(N1*(N2+m))) + 2*log((N1+m-1)/N1) + log((N1+1)/(N1+m)))


ζ = 0.1
m = 10
s(ζ, N) = sum([q[1, i+1] * exp(-ζ^2*Δt*(N-1-i)) for i in 0:N-1])

#####
# Linear iterpolation of s
m = 20
err_linear(ζ, m) = (s.(ζ,1:m) .- ((s(ζ, m) - s(ζ, 1)) / (m-1) .* collect(0:m-1) .+ s(ζ, 1))) 
err_linear(0.1, m)

quadgk(ζ -> maximum(err_linear(ζ, m)) / ζ * sin(ζ*B/rb) * (1- exp(-ζ^2 * constants.Δt̃)), 0, Inf)

plot( (s(m) - s(1)) / (m-1) .* collect(0:m-1) .+ s(1))
plot!(s.(1:m))

bound(block.N[1], block.N[2], m) / 5



#####
# Quadratic interpolation of s
m1 = 4
m2 = 2

x3 = Nt
x2 = x3 - m2
x1 = x3 - m1

xx = x1:x3

s_quad = @. s(x1) * (xx - x2) * (xx - x3) / (x1-x2) / (x1-x3) + s(x2) * (xx - x1) * (xx - x3) / (x2-x1) / (x2-x3) + s(x3) * (xx - x1) * (xx - x2) / (x3-x1) / (x3-x2)

(s.(xx) - s_quad) #/ (2 * π^2 *kg)

plot(s.(xx))
plot!(s_quad)




