using BenchmarkTools
using FiniteLineSource

#export SUITE

SUITE = BenchmarkGroup()

SUITE["non-historical"] = BenchmarkGroup()

function compute_integral_throught_history()
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]

    α, kg = 10^-6, 3.
    r, rb = 0.5, 0.1
    Δt = 3600.
    Δt̃ = α*Δt/rb^2
    C = 1 / (2 * π^2 * r * kg) 

    I = zeros(length(q))
    x, fx, v = precompute_parameters(r/rb)
    compute_integral_throught_history!(I, q, fx, v, C, x, r, kg, Δt̃)
end

function compute_convolution()
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]

    FiniteLineSource.convolve_step(q, Δt = 3600, r = 0.5)
end

function compute_historical_convolution()
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]

    for i in 1:length(q)
        @views FiniteLineSource.convolve_step(q[1:i], Δt = 3600, r = 0.5)
    end
end


SUITE["non-historical"]["compute_integral_throught_history"] = @benchmarkable compute_integral_throught_history()
SUITE["non-historical"]["compute_convolution"] = @benchmarkable compute_convolution()
SUITE["non-historical"]["compute_historical_convolution"] = @benchmarkable compute_historical_convolution()
