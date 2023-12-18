using FiniteLineSource
using Test
using SpecialFunctions

# include("Aqua.jl")

const ϵ = 5*10^-14
# @testset "segment_to_segment_step_response" begin
#     source_segment = BoreholeSegment(0, 0, 1, 5, 0.5)
#     receptor_segment = BoreholeSegment(0, 1, 1, 5, 0.5)
#     @test segment_to_segment_step_response(0, source=source_segment, receptor=receptor_segment)        ≈ 0        atol = 5*10^-20
#     @test segment_to_segment_step_response(100, source=source_segment, receptor=receptor_segment)      ≈ 0        atol = 5*10^-20
#     @test segment_to_segment_step_response(1000, source=source_segment, receptor=receptor_segment)     ≈ 0        atol = 5*10^-10
#     @test segment_to_segment_step_response(100000, source=source_segment, receptor=receptor_segment)   ≈ 0.131890 atol = 5*10^-6
#     @test segment_to_segment_step_response(10000000, source=source_segment, receptor=receptor_segment) ≈ 1.074141 atol = 5*10^-6
# end

@testset "Bakhvalov_Vasil'eva_analytical" begin
    short_step_response(n, t, r, α = 10^-6, kg = 3.) = ( erf( r / sqrt(n*4*α*t) ) - erf( r / sqrt((n+1)*4*α*t) ) ) / (4*π*r*kg) 

    @test compute_integral([1], Δt = 3600, r = 1)               ≈ point_step_response(3600, 1)              atol = ϵ
    @test compute_integral([1], Δt = 3600, r = 0.1)             ≈ point_step_response(3600, 0.1)            atol = ϵ
    @test compute_integral([1; zeros(5)], Δt = 3600, r = 0.1)   ≈ short_step_response(5, 3600, 0.1)         atol = ϵ
end

@testset "Bakhvalov_Vasil'eva_convolution" begin
    error(Q; Δt, r) = compute_integral(Q, Δt = Δt, r = r) - FiniteLineSource.convolve_step(Q, Δt = Δt, r = r)   

    @test error([1], Δt = 3600, r = 1)                  ≈ 0   atol = ϵ
    @test error([1], Δt = 3600, r = 0.1)                ≈ 0   atol = ϵ
    @test error([1; zeros(5)], Δt = 3600, r = 0.1)      ≈ 0   atol = ϵ
end

# integral over an interval is equal to integral over sub-intervals (Currently Failing!)
@test compute_integral([1; zeros(5)], Δt = 3600, r = 0.1, a=0., b=10.)  ≈ compute_integral([1; zeros(5)], Δt = 3600, r = 1, a=0.,b=5.) +  compute_integral([1; zeros(5)], Δt = 3600, r = 1, a=5.,b=10.)

####################
## Longer load tests
##.
q = [20 *sin(2π*i/8760) + 5. for i=1:1000]
# @test compute_integral(q, Δt = 3600., r = 0.1, n = 300) - FiniteLineSource.convolve_step(q, Δt = 3600, r = 0.1) ≈ 0   atol = 5*10^-14

α, kg = 10^-6, 3.
r, rb = 0.5, 0.1
Δt = 3600.
Δt̃ = α*Δt/rb^2

#########################
# # Single interval
#########################

# compute discretization and frequency parameters
dp = discretization_parameters(0.,10.,100)
fp = frequency_parameters(dp,r/rb)
fp = frequency_parameters_2(dp,r/rb)

##.
# update function f over time. The update is divided into two pieces one that happens before computing the integral and one that happens after
fx = zeros(dp.n + 1)
n = 1000
for q in q[1:n-1]  
    FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q)
    compute_integral(q, fx, dp, fp, r, kg)
    FiniteLineSource.fevolve_2!(fx, dp.x, Δt̃, q)
end

FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q[n])
T_laststep = compute_integral(q[n], fx, dp, fp, r, kg)

##.
FiniteLineSource.convolve_step(q[1:n], Δt = 3600, r = 0.5)
compute_integral(q[1:n], Δt = 3600., r = 0.5, n = dp.n) 
# T - FiniteLineSource.convolve_step(q[1:n], Δt = 3600, r = 0.5)
##.

#########################
# # Discretized Intervals
#########################
r = 0.5 

# dv= [0.,0.06,1.,3.,10.] # Discretization intervals limits
# nv = [60,50,40,50]      # Number of discretization points within a limit

# dv= [0.,10.] # Discretization intervals limits
# nv = [70]      # Number of discretization points

dv= [0.,0.1,10.] # Discretization intervals limits
nv = [30,70]      # Number of discretization points
# dv= [0.,0.1,10.] # Discretization intervals limits
# nv = [40]      # Number of discretization points
Nd = length(nv)         # Number of discretization intervals

dps = [discretization_parameters(a,b,n) for (a,b,n) in zip(dv[1:end-1],dv[2:end],nv)] # Discretization parameters for each interval
fps = [frequency_parameters(dp,r/rb) for dp in dps] # Frequency parameters for each interval

##.
fxs = [zeros(dp.n+1) for dp in dps]  # time dependent function for each interval

n = 10
for q in q[1:n-1]  
    for (fx,dp) in zip(fxs,dps)
        FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q)
        FiniteLineSource.fevolve_2!(fx, dp.x, Δt̃, q)
    end
end

for (fx,dp) in zip(fxs,dps)
    FiniteLineSource.fevolve_1!(fx, dp.x, Δt̃, q[n])
end

T = sum(compute_integral(q[n], fx, dp, fp, r, kg) for (fx,dp,fp) in zip(fxs,dps,fps))
# [compute_integral(q[n], fx, dp, fp, r, kg) for (fx,dp,fp) in zip(fxs,dps,fps)]
##.
# function compute_integral_discretized(q, fxs, dps, fps, r, kg)
#     sum(compute_integral(q[n], fx, dp, fp, r, kg) for (fx,dp,fp) in zip(fxs,dps,fps))
# end

# @time compute_integral_discretized(q[n], fxs, dps, fps, r, kg)
FiniteLineSource.convolve_step(q[1:n], Δt = 3600, r = r)
compute_integral(q[1:n], Δt = 3600., r = r, n = 300) 
