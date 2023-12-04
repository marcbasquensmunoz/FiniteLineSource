using LegendrePolynomials
using SpecialFunctions
using FastGaussQuadrature

# Coefficients cₖ of the approximation of f(x) exp(i ω x) from -1 to 1
# in the form sum_(k=0)^n cₖ J_(k+1/2)(ω) / sqrt(ω)
function integral_coefficients(f::Function, n::Integer)
    x, w = gausslegendre(n+1)
    return [(im)^k * (2k + 1) * sqrt(π/2) * sum(w .* Pl.(x, k) .* f.(x)) for k in 0:n]
end

# Integrates 1/sqrt((a(x+1) + b)^2 - 1) exp(i ω x) from -1 to 1
function auxiliary_integral(a, b)
    f(x) = 1 / sqrt((a*(x+1) + b)^2 - 1)
    c = integral_coefficients(f, n)
    ω(ζ) = ζ*r*a/sqrt(κ)
    return ζ -> 1 / sqrt(ω(ζ)) * sum([c[k+1] * besselj(k+1/2, ω(ζ)) for k in 0:n])
end

# Integrates a sin(ζr/sqrt(κ) (a(x+1) + b))/sqrt((a(x+1) + b)^2 - 1) from -1 to 1
function integral(a, b)
    aux = auxiliary_integral(a, b)
    return ζ -> a * (cos(ζ * r / sqrt(κ) * (a+b)) * imag(aux(ζ)) + sin(ζ * r / sqrt(κ) * (a+b)) * real(aux(ζ))) 
end

# Computes the integral needed in the advancing scheme
function integral(z)
    Δ1 = sqrt(1 + (z-D)^2/r^2) - 1
    Δ2 = sqrt(1 + (z-D-H)^2/r^2) - 1
    Δ3 = Δ2 - Δ1 

    if z > D && z < D + H
        return ζ -> integral(Δ1/2, 1)(ζ) + integral(Δ2/2, 1)(ζ)
    else
        if z <= D 
            return ζ -> integral(Δ3/2, Δ1 + 1)(ζ)
        else
            return ζ -> -integral(-Δ3/2, -Δ1 - 1)(ζ)
        end
    end
end

q = [0, 1, 1, 1, 0, 0, 0, 0]
dt = 1
n = 1000
κ = 10^-6
r = 1
D = 1
H = 50

res = integral(0)(1)

G(ζ) = 0
for t in eachindex(q)
    G(ζ) = G(ζ)*exp(-ζ^2*dt) + q[t] / (2 * π^2 * κ) * (1 - exp(-ζ^2*dt)) / ζ * integral(z)(ζ)
end

using Plots
res = integral(1)
ζ = 0:1:50
fz = [res(i) for i in ζ]
plot(ζ, fz)