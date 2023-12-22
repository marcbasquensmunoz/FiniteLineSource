using FiniteLineSource
using Plots

function compute_max_error(r, n)
    Δt = 3600.

    α, rb, kg = 10^-6, 0.1, 3.
    Δt̃ = α*Δt/rb^2

    H = 8760
    D = 24
    t = 1:20H
    q = [20*sin(2*π/H*i) + 5*sin(2*π/D*i) + 5 for i in t]
    C = 1 / (2 * π^2 * r * kg) 

    I = zeros(length(q))
    x, fx, v = precompute_parameters(r/rb, [0., 0.05, 0.1, 1, 7], Int64.(n .* ones(4)))
    compute_integral_throught_history!(I, q, fx, v, C, x, r, kg, Δt̃)
    conv = FiniteLineSource.convolve_step(q, Δt = Δt, r = r)

    error = @. abs(conv - I)
    return findmax(error)
end


R = [0.1, 1, 10, 100]
N = [20, 40, 60, 80, 100]

for n in N
    for r in R
        println("Max error for r=$r, n=$n: $(compute_max_error(r, n))")
    end
end