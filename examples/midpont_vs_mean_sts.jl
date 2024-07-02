using FiniteLineSource
using FiniteLineSource: step_response
using Plots 
using QuadGK

########

function compare_approximations(;Ds, Hs, Dt, Ht)
    h(z) = step_response(t, SegmentToPoint(D=Ds, H=Hs, r=σ, z=z), Constants())

    z = Dt:Ht/500:Dt+Ht

    midpoint_value = h(Dt+Ht/2)
    mean_value = quadgk(z -> h(z), Dt, Dt+Ht)[1]/Ht
    error = round(abs(midpoint_value - mean_value) / max(midpoint_value, mean_value) * 100, sigdigits=3)

    @show midpoint_value
    @show mean_value
    @show "$error %"
    p1 = Plots.plot(z, h.(z), title="Step response in target line", legend=false)

    params = LevelSetParams(D1=Ds, H1=Hs, D2=Dt, H2=Ht, σ=σ, rb=rb)
    h_mean_sts, r_mean = mean_sts_evaluation(params)

    params_mid = MidPointParams(D=Ds, H=Hs, z=Dt+Ht/2, σ=σ, rb=rb)
    h_midpoint, r_mid = midpoint_evaluation(params_mid)

    diff(r) = h_mean_sts(r) - h_midpoint(r)
    predicted_error = abs(quadgk(r -> diff(r)*point_step_response(t, r*rb, α, kg), r_mean[1], r_mean[end])[1])
    @show predicted_error   
    
    p2 = Plots.plot(r_mean, point_step_response.(t, r_mean*rb, α, kg), title="PS response", legend=false)
    p3 = Plots.plot(r_mean, diff.(r_mean), title="Difference of h(r)", legend=false)

    p4 = Plots.plot(r_mid, h_midpoint.(r_mid), title="Midpoint h", legend=false)
    p5 = Plots.plot(r_mean, h_mean_sts.(r_mean), title="Mean h", legend=false)
    
    Plots.plot(p1, p2, p3, p4, p5)
end

σ = 0.1
t = 3600*24*30*12.

rb = 0.1
α = 10^-6
kg = 3.

### Equal depth and lentgh case 
compare_approximations(Ds=0., Hs=10., Dt=0., Ht=10.)
compare_approximations(Ds=0., Hs=100., Dt=0., Ht=100.)
compare_approximations(Ds=0., Hs=1000., Dt=0., Ht=1000.)

### Target below source, equal length
compare_approximations(Ds=0., Hs=20., Dt=20., Ht=20.)
compare_approximations(Ds=0., Hs=20., Dt=50., Ht=20.)
compare_approximations(Ds=0., Hs=20., Dt=100., Ht=20.)

### Target and source overlap
compare_approximations(Ds=0., Hs=50., Dt=25., Ht=50.)
