using BenchmarkTools
#using FiniteLineSource

#export SUITE

SUITE = BenchmarkGroup()

SUITE["utf8"] = BenchmarkGroup(["string", "unicode"])

teststr = "this is a test"

SUITE["utf8"]["replace"] = @benchmarkable replace($teststr, "a" => "b")
SUITE["utf8"]["join"] = @benchmarkable join($teststr, $teststr)
SUITE["utf8"]["plots"] = BenchmarkGroup()
