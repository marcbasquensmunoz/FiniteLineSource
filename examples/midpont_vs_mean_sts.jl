using FiniteLineSource
using Plots 
using QuadGK

########

function compare_approximations_H(;H, r, Δz)
    compare_approximations(Ds=0., Hs=H, Dt=Δz, Ht=H, σ=r)
end

function compare_approximations(;Ds, Hs, Dt, Ht, σ)
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
    
    return error
end

σ = 100.
t = 3600*24*30*12.

rb = 0.1
α = 10^-6
kg = 3.

### Equal depth and lentgh case 
compare_approximations(Ds=0., Hs=10., Dt=0., Ht=10., σ=σ)
compare_approximations(Ds=0., Hs=100., Dt=0., Ht=100., σ=σ)
compare_approximations(Ds=0., Hs=1000., Dt=0., Ht=1000., σ=σ)

### Target below source, equal length
compare_approximations(Ds=0., Hs=20., Dt=20., Ht=20., σ=σ)
compare_approximations(Ds=0., Hs=20., Dt=50., Ht=20., σ=σ)
compare_approximations(Ds=0., Hs=20., Dt=100., Ht=20., σ=σ)

### Target and source overlap
compare_approximations(Ds=0., Hs=50., Dt=25., Ht=50., σ=σ)

compare_approximations_H(H=10., r=σ, Δz=0.)


r = 0.1:0.1:5
z = 0:0.1:50
pr(r) = compare_approximations_H.(H=10., r=r, Δz=0.)
pz(z) = compare_approximations_H.(H=100., r=0.1, Δz=z)

Plots.plot(z, pz.(z))




function diff_non(;t, σ, Δz, H, α)
    diff(t=t/H^2, σ=σ/H, Δz=Δz/H, H=1., α=α)
end

function diff(;t, σ, Δz, H, α)
    h(z) = step_response(t, SegmentToPoint(D=0., H=H, r=σ, z=z), Constants(α=α))
    midpoint_value = h(Δz+H/2)
    mean_value = quadgk(z -> h(z), Δz, Δz+H)[1] / H
    error = round(abs(midpoint_value - mean_value) / max(midpoint_value, mean_value) * 100, sigdigits=3)
    return error
end


α = 10^-6
H = 10.

t = 3600*24*30*12.
σ = 1.
Δz = 0.

diff(t=t, σ=σ, Δz=Δz, H=H, α=α)
diff_non(t=t, σ=σ, Δz=Δz, H=H, α=α)


function diff_nd(;t̃, σ̃, Δz̃, α)
    h(z) = step_response(t̃, SegmentToPoint(D=0., H=1., r=σ̃, z=z), Constants(α=α))
    midpoint_value = h(Δz̃ + 1/2)
    mean_value = quadgk(z -> h(z), Δz̃, Δz̃+1)[1]
    error = midpoint_value == 0. && mean_value == 0. ? 0. : round(abs(midpoint_value - mean_value) / max(midpoint_value, mean_value) * 100, sigdigits=3)
    return error
end

diff_nd(t̃=3600*24*30*12/100, σ̃=1/10, Δz̃=0., α=10^-6)

σ = 0.1:0.1:40
z = 0:0.1:40
surf = [diff_nd(t̃=300000., σ̃=σt, Δz̃=zt, α=10^-6) for σt in σ, zt in z]

WGLMakie.surface(σ, z, surf)

t = 3600/100:3600/100:3600*24*30*12/100
ft(t) = diff_nd(t̃=t, σ̃=1., Δz̃=0., α=10^-6)
Plots.plot(t, ft.(t))



t = 3600*24*30*12.
α = 10^-6
kg = 3.
σ = 0.1
Hs = 10.
Ht = 10.
params = MeanSegToSegEvParams(D1=0., H1=Hs, D2=15., H2=Ht, σ=σ)

lsl, r_mean = FiniteLineSource.level_set_length(params)
Plots.plot(r_mean, lsl.(r_mean))

h_mean_sts, r_mean = mean_sts_evaluation(params)
Plots.plot(r_mean, h_mean_sts.(r_mean))

quadgk(r -> h_mean_sts(r) * point_step_response(t, r, α, kg), r_mean[1], r_mean[end])
hcubature(z -> point_step_response(t, sqrt(σ^2 + (z[1]-z[2])^2), α, kg) / Ht, [0., 0.], [Hs, Ht])
quadgk(z -> step_response(t, SegmentToPoint(D=0., H=Hs, r=σ, z=z), Constants())  / Ht, 0., Ht)



function get_int(f, rs, rf, n)
    xt, w = gausslegendre(n+1)
    x = @. xt * (rf-rs)/2 + (rf+rs)/2
    (rf-rs)/2 * sum(@. f(x) * w)
end


f(r) = h_mean_sts(r) * point_step_response(t, r, α, kg) - erfc(σ/sqrt(4*t*α)) / (2π*kg*sqrt(2*σ*(r-σ)))
An = get_int(f, r_mean[1], r_mean[end], 50)
Di = sqrt(r_mean[end] - r_mean[1]) * erfc(r_mean[1]/sqrt(4α*t)) / (kg*π*sqrt(2*r_mean[1]))
An + Di

quadgk(r -> h_mean_sts(r) * point_step_response(t, r, α, kg), r_mean[1], r_mean[end])
