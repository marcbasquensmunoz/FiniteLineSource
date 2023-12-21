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

# Computes F in the next time step from the previous one
function evolve!(fx, x, Δt, q)
    @. fx = fx * exp(-x^2*Δt) + q * (1 - exp(-x^2*Δt)) / x
end

# I divided the computation of the integrand function evolution into two pieces. 
function fevolve_1!(fx, x, Ct, q)
    @. fx = Ct * (fx - q/x)
end
function fevolve_2!(fx, x, q)
    @. fx = fx  + q / x
end


# Computes the integral of F from 0 to Inf at time t = t * length(Q) using the Bakhvalov and Vasileva method
function compute_integral(Q; n=100, r, Δt, α = 10^-6, rb = 0.1, kg = 3., a =0., b=10.)
    # Total simulation time
    t = Δt * length(Q)
    # The heat wave has not reached points at distance r yet (up to double precision)
    if t < r^2 / (12^2 * α)
        return 0.
    end

    # Computation parameters
    Δt̃ = α/rb^2 * Δt   # Non-dimensional time step
    r̃ = r/rb           # Non-dimensional distance
    ω = r̃              # Oscillation frequency of the integrand

    # Global constant 
    C = 1 / (2 * π^2 * r * kg) 
    
    # a = 0.
    # b = 10.
    m = (b-a)/2
    c = (b+a)/2

    xt, w = gausslegendre(n+1)
    x = zeros(n+1)
    @. x = m * xt + c

    fx = zeros(n+1)

    for q in Q
        evolve!(fx, x, Δt̃, q)
    end

    last_load = Q[length(Q)] 
    @. fx = fx - last_load / x
    
    Ω = ω * m
    Iexp = zeros(n+1)
    for k in 0:n
        Iexp += (im)^(k) * (2k+1) * Bessels.besselj(k+1/2, Ω) * Pl.(xt, k)
    end
    Iexp = sum(w .* Iexp .* fx) * sqrt(π/(2Ω))
    Iexp *= m*exp(im*ω*c)

    Ic = π/2*last_load
    return C * (imag(Iexp) + Ic)
end


# Compute parameters of the discretization that only depends on the interval and the number of integration points 
function discretization_parameters(a,b,n)
    xt,w = gausslegendre(n+1)    
    M = [w[s]*Pl.(xt[s],k) for k=0:n,s=1:n+1]
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, M=M, n=n, xt=xt, w=w)
end

function frequency_parameters_2(dp, ω)
    Ω = dp.m * ω
    K =  dp.m * exp(im * ω * dp.c)
    v = [ sum( [(im)^k *(2k+1) * sqrt(π/(2Ω)) * Bessels.besselj(k+1/2, Ω)*Pl.(dp.xt[s],k) for k =0:dp.n]) * dp.w[s] for s=1:dp.n+1]
    return (ω=ω, Ω=Ω, v=v, K=K)
end

# Computes the integral of F from 0 to Inf at the next time step t = t + Δt using the Bakhvalov and Vasileva method given the current value of the function fx
function compute_integral(fx, dp, fp, r, kg) 
    C = 1 / (2 * π^2 * r * kg) 
    # Iexp = sum(fp.v .* fx) * fp.K
    Iexp = dot(fx , fp.v* fp.K) 
    return C * imag(Iexp)
end


compute_integral_slow(q, r, kg) = q / (4 * π * r * kg) 

function compute(q; t = 3600., r = 1., dv = [0., 10.], nv = [100], α = 10^-6, rb = 0.1, kg = 3.)
    Δt̃ = α*Δt/rb^2
    
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(dv[1:end-1],dv[2:end],nv)] # Discretization parameters for each interval
    fps = [frequency_parameters_2(dp, r/rb) for dp in dps] # Frequency parameters for each interval
    
    fxs = [zeros(dp.n+1) for dp in dps]  # time dependent function for each interval
    
    n = length(q)
    for q in q[1:n-1]  
        for (fx,dp) in zip(fxs,dps)
            fevolve_1!(fx, dp.x, Δt̃, q)
            fevolve_2!(fx, dp.x, q)
        end
    end
    
    for (fx,dp) in zip(fxs,dps)
        fevolve_1!(fx, dp.x, Δt̃, q[n])
    end
    
    return sum(compute_integral(fx, dp, fp, r, kg) for (fx,dp,fp) in zip(fxs,dps,fps)) + compute_integral_slow(q[n], r, kg)
end

function compute_series(q; Δt = 3600., r = 1., n=100, α = 10^-6, rb = 0.1, kg = 3., b=10.)
    Δt̃ = α*Δt/rb^2

    dp = discretization_parameters(0.,b,n)
    fp = frequency_parameters_2(dp,r/rb)

    fx = zeros(dp.n + 1)
    qn = length(q)
    series = []

    for q in q[1:qn-1]  
        FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q)
        append!(series, compute_integral(fx, dp, fp, r, kg) + compute_integral_slow(q, r, kg))
        FiniteLineSource.fevolve_2!(fx, dp.x, q)
    end

    FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q[qn])
    T_laststep = compute_integral(fx, dp, fp, r, kg) + compute_integral_slow(q[qn], r, kg) 
    append!(series, T_laststep)
end



# precompute all the parameters necessary for the computation
function precompute_parameters(r, dv = [0., 10.], nv = [100])
    dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(dv[1:end-1],dv[2:end],nv)] # Discretization parameters for each interval
    fps = [frequency_parameters_2(dp, r) for dp in dps] # Frequency parameters for each interval
    fxs = [zeros(dp.n+1) for dp in dps]  # time dependent function for each interval
    
    x  = reduce(vcat, (dp.x for dp in dps))
    fx = reduce(vcat, fxs)
    v  = reduce(vcat, (fp.v * fp.K for fp in fps))
    
    # return dps,fps,fxs
    return x,fx,v
end

# advance f through the whole load history q
function fevolve_throughhistory!(fx, x, Ct, q)
    for q in q
            fevolve_1!(fx, x, Ct, q)
            fevolve_2!(fx, x, q)
    end
    return nothing
end

# advance f one step and compute integral value 
function compute_integral_and_advance_one_step!(fx, q, v, C, x, Ct, r, kg)
    fevolve_1!(fx, x, Ct, q)
    Iv = C * imag(dot(fx , v)) + compute_integral_slow(q, r, kg) #compute_integral(fx, v, C)
    fevolve_2!(fx, x, q)    
    # Iv += compute_integral_slow(q, r, kg)
    return Iv
end 

function compute_integral_throughthistory!(I, q, fx, v, C, x, r, kg, t)
    N = length(I) 
    Ct = @. exp(-x^2*t)
        for k=1:N 
            # @views @inbounds I[k] = compute_integral_and_advance_one_step!(fx, q[k], v, C, x, r, kg, t)

            @views @inbounds fevolve_1!(fx, x, Ct, q[k])
            @views @inbounds I[k] = C * imag(dot(fx , v)) + compute_integral_slow(q[k], r, kg)  #compute_integral(fx, v, C)
            @views @inbounds fevolve_2!(fx, x, q[k])    
        end
    
       return nothing
end

###################
###### Tests ######
###################

###### Check correctness against analytical solutions ######

# Δt = 3600, rb = 0.1, α = 10^-6, kg = 3

# For Q = [1]
# Ie = - erf( r / sqrt(4*α*Δt) ) / (4*π*r*kg)
# I = erfc( r / sqrt(4*α*Δt) ) / (4*π*r*kg)

# r = 1
# Ie = -0.026525823848649224
# I = 1.2355594071082283e-33

# r = 0.1
# Ie = -0.20196952486650624
# I = 0.06328871361998596


# For Q = [1, zeros(n)]
# I = ( erf( r / sqrt(n*4*α*Δt) ) - erf( r / sqrt((n+1)*4*α*Δt) ) ) / (4*π*r*kg) 

# For n = 5, r = 0.1
# I = 0.008558835281496709


###### Check correctness against convolution ######

point_step_response(t, r, α = 10^-6, kg = 3.) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)

function convolve_step(Q; Δt, r, α = 10^-6)
    Q = diff([0; Q])
    t = Δt:Δt:Δt*length(Q)
    response = point_step_response.(t, r, α)
    return conv(Q, response)[1:length(Q)]
end

function compare(Q, Δt, r)
    println(convolve_step(Q, Δt = Δt, r = r))
    println(compute_integral(Q, Δt = Δt, r = r))
end 

compare([1], 3600, 0.1)