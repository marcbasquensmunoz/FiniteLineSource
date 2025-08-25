using FiniteLineSource
import FiniteLineSource: compute_integral_slow, f_evolve_1!, f_evolve_2!, compute_integral_oscillatory, convolve_step, discretization_parameters, frequency_parameters
using Test
using SpecialFunctions

include("Aqua.jl")
const ϵ = 5*10^-14

@testset "point_step_response" begin
    @test point_step_response(3600, 1, 10^-6, 3)       ≈ 0                     atol = ϵ
    @test point_step_response(3600, 0.5, 10^-6, 3)     ≈ 2.0173737265e-10      atol = ϵ
    @test point_step_response(3600, 0.1, 10^-6, 3)     ≈ 0.06328871361998596   atol = ϵ
    @test point_step_response(3600, 0.05, 10^-6, 3)    ≈ 0.29480258983068475   atol = ϵ
    @test point_step_response(3600, 0.01, 10^-6, 3)    ≈ 2.4037320017694954    atol = ϵ

    @test point_step_response(3600, 0.1, 5*10^-6, 3)   ≈ 0.15866725326935396   atol = ϵ
    @test point_step_response(3600, 0.05, 5*10^-6, 3)  ≈ 0.42024724353889786   atol = ϵ
    @test point_step_response(3600, 0.01, 5*10^-6, 3)  ≈ 2.5410870574168927    atol = ϵ

    @test point_step_response(3600*24, 0.1, 10^-6, 3)  ≈ 0.21483109039800657   atol = ϵ
    @test point_step_response(3600*24, 0.05, 10^-6, 3) ≈ 0.4797249950834168    atol = ϵ
    @test point_step_response(3600*24, 0.01, 10^-6, 3) ≈ 2.6016733120703512    atol = ϵ
end

@testset "convolve_step" begin
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    conv = convolve_step(q, Δt = 0.36, r = 1)
    @test conv[length(q)] ≈ 0.0007359029638686509 atol = ϵ
end