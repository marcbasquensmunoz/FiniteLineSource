using FiniteLineSource

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
I = zeros(length(q))

#setup = PointToPoint(r=1.)
setup = SegmentToPoint(D=10., H=100., σ=0.0575, z=60.)
#setup = SegmentToSegment(D1=0., H1=50., D2=0., H2=50., σ=0.1)
#setup = SegmentToSegmentOld(D1=5., H1=10., D2=5., H2=10., σ=100.)

#setup = MovingPointToPoint(x=100., σ=0.1, v=.01)
#setup = MovingSegmentToPoint(x=0.01, y=0., z=10., v=0.0001, D=0., H=20.)

params = Constants(Δt = 3600.)
precomp = @time precompute_parameters(setup, params=params)
@time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
I

C = convolve_step(q, setup; params=params)
err = abs.(C .- I)
maximum(err)
