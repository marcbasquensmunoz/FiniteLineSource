using FastGaussQuadrature
using Bessels
using LegendrePolynomials
using SpecialFunctions
using DSP

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

# Computes the integral of F from 0 to Inf at time t = Δt * length(Q) using the Bakhvalov and Vasil’eva method
function compute_integral(Q; n=100, r, Δt, α = 10^-6, rb = 0.1, kg = 3.)
    # Total simulation time
    t = Δt * length(Q)
    # The heat wave has not reached points at distance r yet (up to double precision)
    if t < r^2 / (12^2 * α)
        return 0
    end

    # Computation parameters
    Δt̃ = α/rb^2 * Δt   # Non-dimensional time step
    r̃ = r/rb           # Non-dimensional distance
    ω = r̃              # Oscillation frequency of the integrand

    # Global constant 
    C = 1 / (2 * π^2 * r * kg) 

    a = 0
    b = 10
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

point_step(t, r, α = 10^-6, kg = 3.) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)

function convolve_step(Q; Δt, r, α = 10^-6)
    Q = diff([0; Q])
    t = Δt:Δt:Δt*length(Q)
    response = point_step.(t, r, α)
    return conv(Q, response)[length(Q)]
end

function compare(Q, Δt, r)
    println(convolve_step(Q, Δt = Δt, r = r))
    println(compute_integral(Q, Δt = Δt, r = r))
end 

compare([1], 3600, 0.1)

