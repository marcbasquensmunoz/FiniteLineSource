using QuadGK
using FastGaussQuadrature
using LegendrePolynomials
using Bessels

function test_L2(n, a, b, σ, ϵ)
    # Compute the n necessary to meet the tolerance ϵ
    f(x) = 1/sqrt(x^2-σ^2)
    m = (b-a)/2
    c = (b+a)/2
    # Compute the Legendre coefficients and build the Legendre approximation
    ck = [(2k+1)/2 * quadgk(x-> f(m*x+c)*Pl(x, k), -1, 1, atol=ϵ/1000)[1] for k in 0:n]
    fn(x) = sum([Pl((x-c)/m, k)*ck[k+1] for k in 0:n])

    # Compute the L2 norm of f-fn
    # We are divding by m because the error bound is for the interval (-1, 1)
    I, _ = quadgk(x -> 1/m *(f(x)-fn(x))^2, a, b, atol=ϵ/1000)
    sqrt(I)
end

function test_bakhalov(n, a, b, σ, ω, ϵ)
    # Compute the n necessary to meet the tolerance ϵ (on the integral I)
    m = (b-a)/2
    c = (b+a)/2

    # Use the Bakhalov method
    x, w = gausslegendre(n+1)
    P = [w[s] * Pl(x[s], k) for k in 0:n, s in 1:n+1]
    X = @. m*x+c
    f = @. 1 / sqrt(X^2-σ^2)
    Bessels.besselj!(X, 1/2:(n+1/2), m*ω)
    @. X = X * imag(exp(im*ω*c) * im^(0:n)) * (2(0:n)+1)

    # Compute the integral with Bakhalov
    Ib = sqrt(m*π/(2ω)) * X' * P * f
    
    # Compute the integral explicitly
    Ig, E = quadgk(x -> sin(ω*x)/sqrt(x^2-σ^2), a, b, atol=ϵ/1000)
    abs(Ib-Ig)
end


N_test = [10, 25, 50, 75, 100, 125, 150, 200, 500]

ω = 10.

vσ = [1., 5., 10., 20., 50., 100., 200., 400., 700., 1000.]
vb = [2., 4., 6., 10., 20., 35., 50., 75., 100., 200., 500., 1000.]
va = [1., 3., 5., 8., 10., 15., 25., 35., 50., 75., 100.]
vϵ = [10. ^k for k in -4:-2:-12]

min_n = zeros(length(vϵ), length(vσ), length(va), length(vb))

for (i, ϵ) in enumerate(vϵ)
    for (j, σ) in enumerate(vσ)
        for (k, A) in enumerate(va)
            for (l, B) in enumerate(vb)
                b = σ * B
                if k == 1 || σ + 2 > σ * A
                    a = min(σ + max(2, (b-σ) * 0.01), b)
                else 
                    a = σ * A
                end
                @info "Testing ϵ=$ϵ, σ=$σ, a=$a, b=$b"
                if a >= b 
                    min_n[i, j, k, l] = N_test[1]
                    continue
                end
                found = false
                for n_test in N_test
                    #if !found && abs(test_bakhalov(n_test, a, b, σ, ω, ϵ)) < ϵ
                    if !found && abs(test_L2(n_test, a, b, σ, ϵ)) < ϵ
                        min_n[i, j, k, l] = n_test
                        found = true
                    end
                end
            end
        end
    end    
end
