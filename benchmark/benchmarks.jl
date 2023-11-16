using BenchmarkTools
#using FiniteLineSource

const SUITE = BenchmarkGroup()
#=SUITE["seg_to_seg"] = @benchmarkable begin
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
end
=#

SUITE["utf8"] = BenchmarkGroup(["string", "unicode"])

teststr = "this is a test"

SUITE["utf8"]["replace"] = @benchmarkable replace($teststr, "a" => "b")
SUITE["utf8"]["join"] = @benchmarkable join($teststr, $teststr)
SUITE["utf8"]["plots"] = BenchmarkGroup()

SUITE["trigonometry"] = BenchmarkGroup(["math", "triangles"])
SUITE["trigonometry"]["circular"] = BenchmarkGroup()
for f in (sin, cos, tan)
    for x in (0.0, pi)
        SUITE["trigonometry"]["circular"][string(f), x] = @benchmarkable ($f)($x)
    end
end

SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
for f in (sin, cos, tan)
    for x in (0.0, pi)
        SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
    end
end