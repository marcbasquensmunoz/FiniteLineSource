using FiniteLineSource
using FiniteLineSource: analytical_test

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
I = zeros(length(q))

#setup = PointToPoint(r=5.)
#setup = SegmentToPoint(D=5., H=55., r=5., z=30.)
setup = SegmentToSegment(D1=5., H1=20., D2=5., H2=15., r=1.)

params = Constants(Δt = 3600.)
prealloc = Preallocation(setup, params) 

precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
@time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
I    

C = convolve_step(q, setup; params=params)
err = abs.(C .- I)
maximum(err)
