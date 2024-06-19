using FiniteLineSource
using FiniteLineSource: point_step_response
using LegendrePolynomials
using QuadGK
using Plots
using FastGaussQuadrature
using Cubature
using Bessels
using Parameters

function full_function(ck, z; a, b)
    n = length(ck)
    f = zeros(length(z))
    x = @. (2z-a-b)/(b-a)
    for i in 0:n-1
        @. f += ck[i+1] * Pl(x, i)
    end
    f
end

function D_matrix(geometry::SegmentToSegment; N, xN = N, rb=rb)
    @unpack D1, D2, H1, H2, σ = geometry
    x, w = gausslegendre(xN+1)
    π/2 * H1/2 .* [(2l+1)/2 * sum([w[i]*w[j] * Pl(x[i], k) * Pl(x[j], l) / r̃(x[i], x[j], geometry, rb) for i in 1:xN+1, j in 1:xN+1]) for l in 0:N, k in 0:N]
end

function r̃(z1, z2, geometry::SegmentToSegment, rb) 
    @unpack D1, D2, H1, H2, σ = geometry
    sqrt(σ^2 + (H1/2*(z1+1) + D1 - H2/2*(z2+1) - D2)^2) / rb
end
function R_matrix(geometry::SegmentToSegment; n, N, m, c, xN = N, rb=rb)
    @unpack D1, D2, H1, H2, σ = geometry

    x, w = gausslegendre(xN+1)
    # M_lka
    M = zeros(N+1, N+1, n+1)

    R(r, a) = imag(im^a * exp(im*r*c)) / r^(3/2) * besselj(a+1/2, r*m)

    for a in 0:n
        for k in 0:N
            for l in 0:N
                M[l+1, k+1, a+1] = (2l+1)/2 * sum([w[i]*w[j] * Pl(x[i], k) * Pl(x[j], l) * R(r̃(x[i], x[j], geometry, rb), a) for i in 1:xN+1, j in 1:xN+1])
            end
        end
    end
    M .* (H1/2)
end

function P_matrix(;n, m)
    x, w = gausslegendre(n+1)
    sqrt(m*π/2) .* [(2a+1) * w[b] * Pl(x[b], a) for a in 0:n, b in 1:n+1]
end

function advance_1!(f, qk, A, B)
    N = length(qk)
    for k in 1:N
        f[:, k] = A .* f[:, k] + B .* qk[k]
    end
end
function advance_2!(f, q, C)
    N = length(q)
    for k in 1:N
        f[:, k] += C .* qk[k]
    end
end
function integral!(Tk, f; qk, RP, D)
    N = size(RP)[1]-1
    n = size(RP)[3]-1
    mult = reshape(reshape(RP, (N+1)^2, n+1)*f, N+1, N+1, n+1)
    for i in 1:N+1
        Tk[i] = (sum([mult[i, j, j] for j in 1:N+1]) + (D*qk)[i]) .* constant
    end
end

########################################
############### REMEMBER ###############
########################################
## z' is the source line -> [D1, D1+H1]
## z  is the target line -> [D2, D2+H2]
########################################


D1 = 0.
H1 = 1.
D2 = 0.
H2 = 1.
σ = 0.1

α = 1e-6
kg = 3.
Δt = 3600.
rb = 0.1

geometry = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=σ)
params = Constants(α=α, kg=kg, Δt=Δt, rb=rb)
n = 100
N = 10

T = 1

a = 0.
b = 10.
m = (b-a)/2
c = (b+a)/2

constant = 1 / (2 * π^2 * kg * rb)

x, w = gausslegendre(n+1)
ζ = @. m*x + c

f = zeros(n+1, N+1)
Tk = zeros(N+1, T)
X, W = gausslegendre(N+1)
Z = @. H1/2*(X+1) + D1
Q = ones(length(Z))#sin.(Z)
qk = legendre_coeffs(Q, X, W)

Δt̃ = α*Δt/rb^2
# Advance F to next time step
A = @. exp(-Δt̃ * ζ^2)
B = @. -exp(-Δt̃ * ζ^2) / ζ
C = @. 1/ζ

RR = @time R_matrix(geometry, n=n, N=N, m=m, c=c, rb=rb, xN=100)
P = P_matrix(n=n, m=m)
RP = reshape(reshape(R, (N+1)^2, n+1)*P, N+1, N+1, n+1)
D = D_matrix(geometry, N=N, xN=100)

for i in 1:T
    advance_1!(f, qk, A, B)
    integral!(@view(Tk[:, i]), f, qk=qk, RP=RP, D=D)
    advance_2!(f, qk, C)
end

res = full_function(Tk[:, 1], Z, a=D2, b=D2+H2)
plot(Z, res)



# Test advance
#=
q = 1.
fx = zeros(length(ζ))
Ct = @. exp(-ζ^2*Δt̃)
for t in 1:T
    @. fx = Ct * (fx - q/ζ)
    @. fx = fx + q / ζ
end
sum(abs.([full_function(f[i, :], Z, a=D1, b=D1+H1)[1] - fx[i] for i in eachindex(fx)]))
=#

# Test integral
#=
setup = SegmentToPoint(D=D1, H=H1, r=σ, z=D2+H2/2)
precomp = precompute_parameters(setup, params=params)
q = [1.]
I = zeros(length(q))
compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
I
precomp.I_c

I-ff[1]
=#

#precomp.I_c .* ones(length(ff)) + ff

#=
# Test D matrix
res = full_function(D*qk, Z, a=D2, b=D2+H2) * constant
plot(Z, res)

Ic(z) = FiniteLineSource.constant_integral(SegmentToPoint(D=D1, H=H1, r=σ, z=z), params=params)
plot!(Z, Ic.(Z))
=#



function test(z)
    setup = SegmentToPoint(D=D1, H=H1, r=σ, z=z)
    precomp = precompute_parameters(setup, params=params)
    q = [1.]
    I = zeros(length(q))
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
    I[1]
end



n = 100
K = 10
nk = 100
T = 1

Tkk = zeros(K+1, T)

X, W = gausslegendre(K+1)
Z = @. H1/2*(X+1) + D1
Q = ones(K+1)
qk = legendre_coeffs(Q, X, W)

precomp = precompute_matrices(geometry, params, n=n, K=K, nk=nk)
compute_coefficients_through_history(Tkk, qk, params, precomp)

zzz = D2:0.01:D2+H2
plot(zzz, full_function(Tkk[:, 1], zzz, a=D2, b=D2+H2))
plot!(zzz, test.(zzz))
