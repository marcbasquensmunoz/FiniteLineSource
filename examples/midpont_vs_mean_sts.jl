using FiniteLineSource
using FiniteLineSource: point_step_response, segment_step_response, segment_mean_step_response
using Plots 
using QuadGK

rb = 0.1
σ = 5.
α = 10^-6
kg = 3.
t = 3600*24*30*12.

## Source line
D = 0.
H = 20.

## Target line
Dt = 20.
Ht = 20.

####

params = LevelSetParams(D1=D, H1=H, D2=Dt, H2=Ht, σ=σ, rb=rb)
h_mean_sts, r_mean = mean_sts_evaluation(params)

## Comparison of level set method with actual double integral
pred_mean = quadgk(r -> h_mean_sts(r)*point_step_response(t, rb*r, α, kg), r_mean[1]+10^-14, r_mean[end])[1]
exact_mean = segment_mean_step_response(t, params.σ, params.D1, params.H1, params.D2, params.H2, α, kg)
pred_mean - exact_mean

Plots.plot(r_mean, h_mean_sts.(r_mean))
Plots.plot(r_mean, point_step_response.(t, rb.*r_mean, α, kg))
####

params_mid = MidPointParams(D=D, H=H, z=Dt+Ht/2, σ=σ, rb=rb)
h_midpoint, r_mid = midpoint_evaluation(params_mid)

## Comparison of midpoint method with actual double integral
pred_mid = quadgk(r -> h_midpoint(r)* point_step_response(t, r*rb, α, kg), r_mid[1], r_mid[end])[1]
exact_mid = segment_step_response(t, params_mid.σ, params_mid.z, params_mid.D, params_mid.H, α, kg)
pred_mid - exact_mid

Plots.plot(r_mid, h_midpoint.(r_mid))
####

## Comparison of two approximations
diff(r) = h_mean_sts(r) - h_midpoint(r)
Plots.plot(r_mean, diff.(r_mean))

quadgk(r -> diff(r) * point_step_response(t, r*rb, α, kg), r_mean[1]+10^-14, r_mean[end])[1]



########

function compare_approximations(;Ds, Hs, Dt, Ht)
    z = Dt:Ht/500:Dt+Ht
    h(z) = segment_step_response(t, σ, z, Ds, Hs, α, kg)

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
