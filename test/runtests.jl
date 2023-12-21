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

@testset "compute_integral_slow" begin
    @test compute_integral_slow(1, 1, 3)   ≈ 0.026525823848649224   atol = ϵ
    @test compute_integral_slow(10, 1, 3)  ≈ 0.26525823848649224    atol = ϵ
    @test compute_integral_slow(100, 1, 3) ≈ 2.6525823848649224     atol = ϵ

    @test compute_integral_slow(1, 5, 3)   ≈ 0.005305164769729845   atol = ϵ
    @test compute_integral_slow(10, 5, 3)  ≈ 0.05305164769729845    atol = ϵ
    @test compute_integral_slow(100, 5, 3) ≈ 0.5305164769729845     atol = ϵ
end

@testset "discretization_parameters" begin
    @testset "0_10_100" begin
        dp = discretization_parameters(0, 10, 100)
        @test dp.m == 5.0        
        @test dp.c == 5.0        
        @test dp.n == 100
        @test length(dp.w) == 101 && length(dp.x) == 101 && length(dp.xt) == 101 
        @test findmin(dp.xt)[1] > -1 && findmax(dp.xt)[1] < 1
        @test sort(dp.xt) == dp.xt
        @test findmin(dp.x)[1] > 0 && findmin(dp.x)[1] < 0.1 && findmax(dp.x)[1] < 10 &&findmax(dp.x)[1] > 9.9
        @test sort(dp.x) == dp.x
        @test size(dp.M) == (101, 101)
    end 

    @testset "2_7_200" begin
        dp = discretization_parameters(2, 7, 200)
        @test dp.m == 2.5        
        @test dp.c == 4.5        
        @test dp.n == 200
        @test length(dp.w) == 201 && length(dp.x) == 201 && length(dp.xt) == 201 
        @test findmin(dp.xt)[1] > -1 && findmax(dp.xt)[1] < 1
        @test sort(dp.xt) == dp.xt
        @test findmin(dp.x)[1] > 2 && findmin(dp.x)[1] < 2.1 && findmax(dp.x)[1] < 7 &&findmax(dp.x)[1] > 6.9
        @test sort(dp.x) == dp.x
        @test size(dp.M) == (201, 201)
    end 
end

@testset "frequency_parameters" begin
    @testset "0_10_100_10" begin
        dp = discretization_parameters(0, 10, 100)
        fp = frequency_parameters(dp, 10.)
        @test fp.ω == 10.
        @test fp.Ω == 50.
        @test length(fp.v) == 101
        @test fp.K ≈ 4.824830142460566 - 1.3118742685196438im  atol = ϵ
    end
end

@testset "precompute_parameters" begin
    @testset "10" begin
        x, fx, v = precompute_parameters(10)
        @test length(fx) == 101
        @test length(x) == 101
        @test length(v) == 101

        @test findmin(fx)[1] == 0 && findmax(fx)[1] == 0 
        @test findmin(x)[1] > 0 &&  findmin(x)[1] < 0.1 && findmax(x)[1] < 10 && findmax(x)[1] > 9.9 
        @test sort(x) == x
    end

    @testset "10, [0., 1., 10.], [100, 100]" begin
        x, fx, v = precompute_parameters(10, [0., 1., 10.], [100, 100])
        @test length(fx) == 202
        @test length(x) == 202
        @test length(v) == 202

        @test findmin(fx)[1] == 0 && findmax(fx)[1] == 0 
        @test @views findmin(x[1:101])[1] > 0 &&  findmin(x[1:101])[1] < 0.1 && findmax(x[1:101])[1] < 1 && findmax(x[1:101])[1] > 0.9
        @test @views findmin(x[102:202])[1] > 1 &&  findmin(x[102:202])[1] < 1.1 && findmax(x[102:202])[1] < 10 && findmax(x[102:202])[1] > 9.9 
        @test sort(x) == x
    end
end

@testset "f_evolve_1!" begin
    @testset "1" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        old_fx = copy(fx)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_1!(fx, x, Ct, 1)
        @test old_x == x
        @test old_fx != fx
    end
end

@testset "f_evolve_2!" begin
    @testset "1" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        old_fx = copy(fx)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_2!(fx, x, 1)
        @test old_x == x
        @test old_fx != fx
    end
end

@testset "compute_integral_oscillatory" begin
    @testset "im .* ones(100) ones(100) 1" begin
        fx = im .* ones(100)
        v = ones(100)
        C = 1
        @test compute_integral_oscillatory(fx, v, C) == -100.
    end

    @testset "im .* ones(100) ones(100) 10" begin
        fx = im .* ones(100)
        v = ones(100)
        C = 10
        @test compute_integral_oscillatory(fx, v, C) == -1000.
    end

    @testset "im .* [0:100] [100:-1:0] 1" begin
        fx = im .* [0:100]
        v = [100:-1:0]
        C = 1
        @test compute_integral_oscillatory(fx, v, C) == -166650.
    end
end

@testset "convolve_step" begin
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    conv = convolve_step(q, Δt = 0.36, r = 1)
    @test conv[length(q)] ≈ 0.0007359029638686509 atol = ϵ
end
