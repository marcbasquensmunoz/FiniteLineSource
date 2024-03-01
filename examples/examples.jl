using FiniteLineSource
using FiniteLineSource: analytical_test

q = ones(1)
I = zeros(length(q))

#setup = PointToPoint(r=5.)
#setup = SegmentToPoint(D=5., H=55., r=5., z=30.)
setup = SegmentToSegment(D1=5., H1=20., D2=10., H2=15., r=5.)

params = Constants(Δt = 3600*24*30.)
prealloc = Preallocation(setup, params) 

@btime precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
@time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
I                                            

res = q[1] * analytical_test(setup, params=params, t=params.Δt)[1]
err = res - I[1]


