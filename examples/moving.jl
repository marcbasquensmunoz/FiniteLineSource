using FiniteLineSource
using QuadGK
using Plots

function u(ζ, r, v, rb, α, kg)
    if ζ < abs(v)/2
        0
    else
        rb^2 / α * α/kg * exp(r*v/2) / ( 2π^2*rb^3*r*ζ ) * sin( r * sqrt( ζ^2 - v^2/4) )
    end
end

### Integration of u to obtain a good numerical result when integrating F
u_int(r, v, rb, α, kg) =  1 / (4π*rb*r*kg)

Δt = 3600.
v = 0.00001
r = 0.1

rb = 0.1
α = 10^-6
kg = 3.

ṽ = v * rb / α
r̃ = r/rb
Δt̃ = Δt*α/rb^2


quadgk(ζ -> (-exp(-ζ^2 * Δt̃)) * u(ζ, r̃, ṽ, rb, α, kg), 0, Inf)[1] + u_int(r̃, ṽ, rb, α, kg) 
quadgk(ζ -> exp(-ζ^2 * Δt̃) * (1-exp(-ζ^2 * Δt̃)) * u(ζ, r̃, ṽ, rb, α, kg), 0, Inf)
quadgk(ζ -> exp(-ζ^2 *2* Δt̃) * (1-exp(-ζ^2 * Δt̃)) * u(ζ, r̃, ṽ, rb, α, kg), 0, Inf)


x = 0:0.1:200
Plots.plot(x, @.  f(x, r̃, ṽ, rb, α, kg))
Plots.plot!(x, @. (1-exp(-x^2 * Δt̃)) * f(x, r̃, ṽ, rb, α, kg))
Plots.plot!(x, @. exp(-ζ^2 * Δt̃) *(1-exp(-x^2 * Δt̃)) * f(x, r̃, ṽ, rb, α, kg))


FiniteLineSource.convolve_step_moving([1, 0, 0]; Δt, r, v, α, kg)
FiniteLineSource.moving_point_step_response(Δt, r, v, α, kg)


u_ps(ζ, r, rb, α, kg) =  rb^2 / α * α/kg / ( 2π^2*rb^3*r*ζ ) * sin(r*ζ)

quadgk(ζ -> (1-exp(-ζ^2 * Δt̃)) * u_ps(ζ, r̃, rb, α, kg) , 0, Inf)
FiniteLineSource.convolve_step([1, 0]; Δt, r, α, kg)


x = 0:0.1:10
@gif for i in 0:00
    Plots.plot(x, @. FiniteLineSource.moving_point_step_response(Δt*i, x, 0.001, α, kg))
end


x = 0:0.1:10
@gif for i in 0:500
    Plots.plot(x, @. FiniteLineSource.point_step_response(Δt*i, x, α, kg))
end



### Test

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
I = zeros(length(q))
Δt = 3600.

segment_limits = [0., 0.01, 0.1, 0.5, 1., 3., 10.]
segment_points = 4 .* [10, 25, 10, 10, 10, 10]

setup = MovingPointToPoint(x=0.01, σ=0.0, v=0.05)
params = Constants(Δt=Δt, segment_limits=segment_limits, segment_points=segment_points)
prealloc = Preallocation(setup, params)  
precomp = precompute_parameters(setup, prealloc=prealloc, params=params)
compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

C = FiniteLineSource.convolve_step(q, setup, params=params)

error = @. C - I 
Plots.plot(1:length(q), error)

