using BenchmarkTools
using FiniteLineSource
using Plots

SUITE = BenchmarkGroup()
SUITE["seg_to_seg"] = @benchmarkable begin
    s1 = BoreholeSegment(0, 0, 0.5, 2, 0.1)
    s2 = BoreholeSegment(0, 1, 1, 5, 0.1)

    t = 0:10^3:10^7
    response(x) = segment_to_segment_step_response(x, source=s1, receptor=s2)
    result = response.(t)
end

SUITE["T_field"] = @benchmarkable begin
    s = BoreholeSegment(0, 0, 0.5, 2, 0.1)
    points = [(i, j) for i in 0:0.1:2, j in 0.1:0.1:2]
    points = T_ls.(last.(points), first.(points), 10, Ref(s))
    heatmap(points)
end

#tune!(suite)
#results = run(suite, verbose = true, seconds = 10)
#"Hello"