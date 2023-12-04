# Algorithm stable when n <= ω
function integral(f::Function, n::Integer, ω)
    α(i) = 2i/(im*ω) 
    x = chebpoints(n, -1, 1)
    c = chebinterp(f.(x), -1, 1)
    a(i) = i < length(c.coefs) ? c.coefs[i+1] : 0
    d = [zeros(n)im; a(n); 0]
    for i in n-1:-1:0
        d[i+1] = a(i) - a(i+2) + d[i+3] - α(i+1) * d[i+2]
    end
    I = 2*sin(ω)/ω * sum([d[i] for i in 1:2:n+1]) - im * 2 * cos(ω)/ω * sum([d[i] for i in 2:2:n+1])
end

# Test function
f(x) = exp(16 *(x-1))
Q(ω) = 2 *exp(-16) * sinh(16 + im*ω)/(16 + im*ω)

#β = 10
#f(x) = exp(2π*im*x*β)
#Q(ω) = 2 * sin(2π*β + ω) / (2π*β + ω)

#f(x) = (1 - x^2)^(3/2)
#Q(ω) = 3π/ω^2 * besselj(2, ω)

ω = 5000

m = 1:10:101

# Test results
error = [abs(integral(f, i, ω) - Q(ω)) for i in m]
plot(m, error)
