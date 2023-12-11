using LegendrePolynomials
using SpecialFunctions
using FastGaussQuadrature
using Bessels

# Coefficients cₖ of the approximation of f(x) exp(i ω x) from -1 to 1
# in the form sum_(k=0)^n cₖ J_(k+1/2)(ω) / sqrt(ω)
function integral_coefficients(f::Function, n::Integer)
    x, w = gausslegendre(n+1)
    k -> (im)^k * (2k + 1) * sqrt(π/2) * sum(w .* Pl.(x, k) .* f.(x))
    #return [(im)^k * (2k + 1) * sqrt(π/2) * sum(w .* Pl.(x, k) .* f.(x)) for k in 0:n]
end

# Integrates 1/sqrt((a(x+1) + b)^2 - 1) exp(i ω x) from -1 to 1
function auxiliary_integral(a, b; real)
    f(x) = 1 / sqrt((a*x + a + b)^2 - 1)
    c = integral_coefficients(f, n)
    ω(ζ) = ζ*r*a/sqrt(κ)
    return ζ -> sum([c(k) * Bessels.besselj(k+1/2, ω(ζ)) for k in (real ? 0 : 1):2:n]) / (real ? 1 : im)
end

# Integrates a sin(ζr/sqrt(κ) (a(x+1) + b))/sqrt((a(x+1) + b)^2 - 1) from -1 to 1
function integral(a, b)
    aux_re = auxiliary_integral(a, b, real=true)
    aux_im = auxiliary_integral(a, b, real=false)
    return ζ -> a * κ^(1/4) / sqrt(ζ*r*a) * (cos(ζ * r / sqrt(κ) * (a+b)) * aux_im(ζ) + sin(ζ * r / sqrt(κ) * (a+b)) * aux_re(ζ)) 
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

q = [0; 1; 1; 1; zeros(4)]
dt = 1
n = 100;
κ = 10^-6;
r = 2;
D = 1;
H = 50;

res = integral(25)(1)

T = zeros(length(q))
g = ζ -> 0
for t in 1:length(q)
    g_old = g
    g = ζ -> g_old(ζ)*exp(-ζ^2*dt) + q[t] / (2 * π^2 * κ) * (1 - exp(-ζ^2*dt)) / ζ * real(integral(z)(ζ))
    T[t] = quadgk(ζ -> g(ζ), 0, Inf)[1]
end
t = 1:length(q)
plot(t, T.(t))



using Plots
res = integral(1)
ζ = 0:1:50
fz = [res(i) for i in ζ]
plot(ζ, fz)

@btime integral(-1)(1)

f(x) = 1 / sqrt((Δ1*(x+1) + 1)^2 - 1)
@btime integral_coefficients(f, 100)(1)

z = 5
Δ1 = sqrt(1 + (z-D)^2/r^2) - 1
Δ2 = sqrt(1 + (z-D-H)^2/r^2) - 1
Δ3 = Δ2 - Δ1 
@btime auxiliary_integral(Δ1/2, 1, real=true)(1)


@btime begin
x, w = gausslegendre(100)
sum(w .* Pl.(x, 1) .* f.(x))
end

n = 100
a = ( sqrt(1 + (z-D)^2/r^2) - 1)/2
b = 1
f(x) = 1 / sqrt((a*(x+1) + b)^2 - 1)
c = integral_coefficients(f, n)
ω(ζ) = ζ*r*(sqrt(1 + (z-D)^2/r^2) - 1)/2/sqrt(κ)
@btime sum([c(k) * Bessels.besselj(k+1/2, ω(1)) for k in 0:2:n])

@btime integral(a, b)(1)