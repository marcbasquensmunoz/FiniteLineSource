using BenchmarkTools
using FiniteLineSource

#export SUITE

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
Δt = 3600.
r = 0.5

SUITE = BenchmarkGroup()

SUITE["non-historical"] = BenchmarkGroup()

function compute_integral_throught_history(q, Δt, r)
    segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
    segment_points = [10, 25, 10, 10, 10, 10]

    params = Constants(Δt=Δt, segment_limits=segment_limits, segment_points=segment_points)
    setup = PointToPoint(r=r)
    prealloc = Preallocation(setup, params) 

    precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
end

function compute_convolution(q, Δt, r)
    FiniteLineSource.convolve_step(q, Δt = Δt, r = r)
end

function compute_historical_convolution(q, Δt, r)
    for i in 1:length(q)
        @views FiniteLineSource.convolve_step(q[1:i], Δt = Δt, r = r)
    end
end

SUITE["non-historical"]["compute_integral_throught_history"] = @benchmarkable compute_integral_throught_history(q, Δt, r)
SUITE["non-historical"]["compute_convolution"] = @benchmarkable compute_convolution(q, Δt, r)
SUITE["non-historical"]["compute_historical_convolution"] = @benchmarkable compute_historical_convolution(q, Δt, r)

tune!(SUITE)
results = run(SUITE)