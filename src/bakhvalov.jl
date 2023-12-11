using Plots
using FastGaussQuadrature
using Bessels
using LegendrePolynomials

# Physical parameters
Δt = 3600    # time step
κ = 10^-6    # thermal diffusivity
r = 1        # distance
rb = 0.1     # borehole radius

# Computation parameters
Δt̃ = κ/rb^2 * Δt   # Non-dimensional time step
r̃ = r/rb           # Non-dimensional distance
ω = r̃              # Oscillation frequency
n = 100            # Number of points to evaluate

# Global constant 
C = 1 / (2 * π^2 * r̃ * rb^3)

function evolve!(fx, x, Δt, q)
    @. fx = fx * exp(-x^2*Δt) + q * (1 - exp(-x^2*Δt)) / x
end

function compute_integral(Q)
    n = 100
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
    Iexp = 0
    for k in 0:n
        Iexp += (im)^(k) * (2k+1) * sqrt(π/(2Ω)) * Bessels.besselj(k+1/2, Ω) * sum(w .* Pl.(xt, k) .* fx)
    end
    Iexp *= m*exp(im*ω*c)

    Ic = 0#π/2*last_load
    return C * (imag(Iexp) + Ic)
end

t = 0:10*8760
Q = [rand([10.,20.,30.,45.,0., -10., -20., -30.])  + 60*sin(1/8760* 2π*i) .- 30 for i in t]  

compute_integral(Q)


###### Tests ######

# For Q = [1], Δt = 3600, κ = 10^-6, r = 1, rb = 0.1 
# Ie = -1/(4*π*r*rb^2) * erf( r / (2*sqrt(κ*Δt)) ) = -7.957747154594767
# I = 1/(4*π*r*rb^2) * erfc( r / (2*sqrt(κ*Δt)) ) = 3.7066782213246847e-31  # Not enough precision


point_impulse(t, r, κ=10^-6) = 1/(4*π*κ*t)^(3/2) * exp(-r^2/(4κ*t))
point_step(t, r, κ=10^-6) = erfc(r/(2*sqrt(t*κ))) / (4*π*r*κ)

function simulate(Q)
    Q = [0; Q]
    t = dt * length(Q)

    result = 0
    for i in 1:length(Q)-1
        result += (Q[i+1] - Q[i]) * point_step(t - i*dt, R)
    end
    return result
    #=
    dt = 3600
    t = dt:dt:dt*length(Q)
    response = point_impulse.(t, R)
    return sum(conv(Q, response))
    =#
end
