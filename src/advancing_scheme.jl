using FastChebInterp
using Polynomials
using Plots
using QuadGK


# Compute the oscillatory integral of f(x) e^(iωx) from -1 to 1
# https://doi.org/10.1016/j.cam.2016.10.034

δ(N, i) = [(j == i ? 0.5 : 1) for j in 0:N]
λ(i) = i == 0 ? 1/2 : (-1)^i
α(i) = 2i/(im*ω) 

f(x) = exp(16 *(x-1))
Q(ω) = 2 *exp(-16) * sinh(16 + im*ω)/(16 + im*ω)
#f(x) = 1/(x^2 + 1/16)
ω = 50
n = 100
#ϵ = 10^-16


x = chebpoints(n, -1, 1)
c = chebinterp(f.(x), -1, 1)
a(i) = i < length(c.coefs) ? c.coefs[i+1] : 0


# Algorithm stable when n <= ω
d = [zeros(n)im; a(n); 0]
for i in n-1:-1:0
    d[i+1] = a(i) - a(i+2) + d[i+3] - α(i+1) * d[i+2]
end
d[1] = d[1] - real(d[1]/2)

# Check solutions
y = -1:0.01:1
p_n = ChebyshevT(a.(0:n))
#=
plot(y, p_n.(y), label = "pₙ")
plot!(y, f.(y), label = "f")
plot!(y, f.(y) - p_n.(y), label = "Error")
error = quadgk(x -> f(x) - p_n(x), -1, 1)[1]
=#
error_squared = println("Error squared of p_n $(sum((f.(y) - p_n.(y)).^2))")

Aϕ(x) = ω / (ω - 16*im) * exp(-32 - im*(1+x)*ω) * (-1 + exp((1+x)*(16+ im*ω)))
ϕ = ChebyshevT(d)
errorϕ = Aϕ.(y) - ϕ.(y)
dϕ = derivative(ϕ)
errorϕ = dϕ + im * ω * ϕ - im * ω * p_n
println("Error squared of ϕ $(sum(abs.(errorϕ.(y)).^2))")

plot(y, abs.(errorϕ), label = "Accuracy of ϕ")


I = exp(im*ω)/(im*ω)*(sum(d)) - exp(-im*ω)/(im*ω)*sum(λ.(0:n+1) .* d)
I = 2*sin(ω)/ω * sum([d[i] for i in 1:2:n+1]) - im * 2 * cos(ω)/ω * sum([d[i] for i in 2:2:n+1])

# Legendre Polynomials method

x, w = gausslegendre(n+1)
c = [(2k + 1)/2 * sum(w .* Pl.(x, k) .* f.(x)) for k in 0:n]
I = sum(c .* [ (im)^k * sqrt(2π /ω) * besselj(k+1/2, ω) for k in 0:n])
println("Error of the integral $(abs(I - Q(ω)))")

# Algorithm for n > ω
M = floor(ω)
μ = [1/2; im/ω - 1; zeros(M-1)im]
for i in 2:M 
    μ[i+1] = λ(i) - α(i) * μ[i] + μ[i-1]
end
τ = (λ(M+1) + μ[M]) / μ[M+1]

z = [[a(j) - a(j+2) for j in 0:M-1]; 0im; 0im]
for j in 0:M-1
    z[M+1] += μ[j+1]*z[j+1]
end
z[M+1] *= -1/μ[M+1]
z[M+2] = a(M) - a(M+2) - z[M+1]
h = [1/μ[M+1], -1/μ[M+1]] # Starts with h_M => N -> N - M + 1

v = 0
S1 = 0
S2 = 0
ψ = 0
u = [α(M+1) - τ] # Starts with u_(M+1) => N -> N - (M+1) + 1
d_last = (z[M+2] - ψ * h[2]) / u[1]

N = M+1
while abs(d_last) > ϵ
    N +=1 
    l = 1/u[N - (M+1)]
    push!(u, α(N) + l)
    push!(h,  -l * z[N])
    push!(z, a(N-1) - a(N+1) - l * z[N])
    v = (λ(N) + v) / u[N - M]
    S1 = S1 + v * z[N+1]
    S2 = S2 + v * h[N+1 - M]
    ψ = S1 / (1 + S2)
    d_last = (z[N+1] - ψ * h[N+1 - M]) / u[N - M]
end

d = [zeros(N); d_last]

for i in N-1:-1:0 
    if i >= M + 1
        d[i+1] = (z[i+1] - ψ * h[i+1 - M] + d[i+2]) / u[i - M]
    elseif i == M
        d[i+1] = z[i+1] - ψ * h[i+1 - M] - τ * d[i+2]
    else
        d[i+1] = z[i+1] - α(i+1) * d[i+2] + d[i+3]
    end
end

I = exp(im*ω) / (im*ω) * (sum(d))
Q(ω)
Q_N = ChebyshevT(real.(d)) - ChebyshevT([real(d[1])]) / 2
