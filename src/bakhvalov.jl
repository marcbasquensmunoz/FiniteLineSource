using Bessels
using DSP
using FastGaussQuadrature
using LegendrePolynomials
using SpecialFunctions
using LinearAlgebra

# Physical parameters
# Δt  - time step
# r   - distance to the borehole
# rb  - borehole radius
# α   - thermal diffusivity
# kg  - thermal conductivity

function f_evolve_1!(fx, x, Ct, q)
    @. fx = Ct * (fx - q/x)
end

function f_evolve_2!(fx, x, q)
    @. fx = fx + q / x
end

compute_integral_slow(q, r, kg) = q / (4 * π * r * kg) 
compute_integral_oscillatory(fx, v, C) = C * imag(dot(fx , v))

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end

function frequency_parameters(dp, ω)
    Ω = dp.m * ω
    K = dp.m * exp(im * ω * dp.c)
    v = [ sum( [(im)^k *(2k+1) * sqrt(π/(2Ω)) * Bessels.besselj(k+1/2, Ω)*Pl.(dp.xt[s],k) for k =0:dp.n]) * dp.w[s] for s=1:dp.n+1]
    return (ω=ω, Ω=Ω, v=v, K=K)
end

# Precompute all the parameters necessary for the computation
function precompute_parameters(r, dv = [0., 10.], nv = [100])
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(dv[1:end-1],dv[2:end],nv)] # Discretization parameters for each interval
    fps = [frequency_parameters(dp, r) for dp in dps] # Frequency parameters for each interval
    fxs = [zeros(dp.n+1) for dp in dps]  # time dependent function for each interval
    
    x  = reduce(vcat, (dp.x for dp in dps))
    fx = reduce(vcat, fxs)
    v  = reduce(vcat, (fp.v * fp.K for fp in fps))
    
    return x, fx, v
end

# Advance f through the whole load history q
function fevolve_through_history!(fx, x, Ct, q)
    for q in q
        f_evolve_1!(fx, x, Ct, q)
        f_evolve_2!(fx, x, q)
    end
    return nothing
end

# Advance f one step and compute integral value 
function compute_integral_and_advance_one_step!(fx, q, v, C, x, Ct, r, kg)
    f_evolve_1!(fx, x, Ct, q)
    Iv = compute_integral_oscillatory(fx, v, C) + compute_integral_slow(q, r, kg) 
    f_evolve_2!(fx, x, q)    
    return Iv
end 

# Advance f through the whole load history q and record each integral value in I
function compute_integral_throught_history!(I, q, fx, v, C, x, r, kg, t)
    Ct = @. exp(-x^2*t)
    for k in eachindex(I)
        @views @inbounds f_evolve_1!(fx, x, Ct, q[k])
        @views @inbounds I[k] = compute_integral_oscillatory(fx, v, C) + compute_integral_slow(q[k], r, kg)
        @views @inbounds f_evolve_2!(fx, x, q[k])    
    end
    
    return nothing
end

# Step response of the point source
point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)

# Solution form the discrete convolution for the load q
function convolve_step(q; Δt, r, α = 10^-6, kg = 3.)
    q = diff([0; q])
    t = Δt:Δt:Δt*length(q)
    response = point_step_response.(t, r, α, kg)
    return conv(q, response)[1:length(q)]
end
