using BenchmarkTools
using FiniteLineSource

SUITE = BenchmarkGroup()

function seg_to_seg()
    s1 = BoreholeSegment(0, 0, 0.5, 2, 0.1)
    s2 = BoreholeSegment(0, 1, 1, 5, 0.1)

    t = 0:10^3:10^7
    response(x) = segment_to_segment_step_response(x, source=s1, receptor=s2)
    result = response.(t)
end

function T_field()
    s = BoreholeSegment(0, 0, 0.5, 2, 0.1)
    points = [(i, j) for i in 0:0.1:2, j in 0.1:0.1:2]
    points = T_ls.(last.(points), first.(points), 10, Ref(s))
end

SUITE["fls"]["seg2seg"] = @benchmarkable seg_to_seg()
SUITE["fls"]["T_field"] = @benchmarkable T_field()