using BenchmarkTools
using FiniteLineSource
using FiniteLineSource: PointSource, LineSource

#export SUITE

Δt = 3600.
Nt = 8760*20

ϵ = 1e-6

rb = .1
α  = 1e-6
kg = 3.
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:Nt]

point_to_point = PointToPoint(r = 1.)
line_to_point = SegmentToPoint(D = 0., H = 100., z = 50., σ = 1.)
line_to_line = SegmentToSegment(D1 = 0., H1 = 100., D2 = 0., H2 = 100., σ = 1.)

point_sources = [PointSource(0., 0., 0., rb), PointSource(1., 0., 0., rb)]
line_sources = [LineSource(0., 0., 0., 100., rb), LineSource(1., 0., 0., 100., rb)]

containers = FiniteLineSource.AsymptoticContainers(10)

function compute_integral_throught_history(q, setup, constants, ϵ)
    I = zeros(length(q))
    precomp = precompute_parameters(setup, params=constants)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=constants)
end

function compute_convolution(q, setup, constants, ϵ)    
    convolve_step(q, setup, params=constants, ϵ=ϵ)
end

function blocks_implementation(q, setup, constants, ϵ, positions, containers) 
    Nt = length(q)
    block = prepare_containers(setup, positions, ϵ, Nt, constants, containers, compute_first_block=false);
    Ib = zeros(length(positions), Nt)
    evolve!(Ib, hcat(q, q)', block)
end


SUITE = BenchmarkGroup()

SUITE["point_to_point"] = BenchmarkGroup()
SUITE["point_to_point"]["original-non-history"] = @benchmarkable compute_integral_throught_history(q, point_to_point, constants, ϵ)
SUITE["point_to_point"]["non-history"] = @benchmarkable blocks_implementation(q, point_to_point, constants, ϵ, point_sources, containers)
SUITE["point_to_point"]["convolution"] = @benchmarkable compute_convolution(q, point_to_point, constants, ϵ)

SUITE["line_to_point"] = BenchmarkGroup()
SUITE["line_to_point"]["original-non-history"] = @benchmarkable compute_integral_throught_history(q, line_to_point, constants, ϵ)
SUITE["line_to_point"]["non-history"] = @benchmarkable blocks_implementation(q, line_to_point, constants, ϵ, line_sources, containers)
SUITE["line_to_point"]["convolution"] = @benchmarkable compute_convolution(q, line_to_point, constants, ϵ)

SUITE["line_to_line"] = BenchmarkGroup()
SUITE["line_to_line"]["original-non-history"] = @benchmarkable compute_integral_throught_history(q, line_to_line, constants, ϵ)
SUITE["line_to_line"]["non-history"] = @benchmarkable blocks_implementation(q, line_to_line, constants, ϵ, line_sources, containers)
SUITE["line_to_line"]["convolution"] = @benchmarkable compute_convolution(q, line_to_line, constants, ϵ)

tune!(SUITE)
results = run(SUITE)
