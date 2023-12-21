using FiniteLineSource
import FiniteLineSource: compute_integral_slow, f_evolve_1!, f_evolve_2!, compute_integral_oscillatory, convolve_step, discretization_parameters, frequency_parameters
using Test
using SpecialFunctions

include("Aqua.jl")
include("test_values.jl")

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
        @test dp.w == w_100
        @test dp.xt == xt_100
        @test dp.x == x_100_0_10
        @test size(dp.M) == (101, 101)
    end 

    @testset "2_7_200" begin
        dp = discretization_parameters(2, 7, 200)
        @test dp.m == 2.5        
        @test dp.c == 4.5        
        @test dp.n == 200
        @test dp.w == w_200
        @test dp.xt == xt_200
        @test dp.x == x_200_2_7
        @test size(dp.M) == (201, 201)
    end 
end

@testset "frequency_parameters" begin
    @testset "0_10_100_10" begin
        dp = discretization_parameters(0, 10, 100)
        fp = frequency_parameters(dp, 10.)
        @test fp.ω == 10.
        @test fp.Ω == 50.
        @test fp.v ≈ v_0_10_100_10 atol = ϵ
        @test fp.K ≈ 4.824830142460566 - 1.3118742685196438im  atol = ϵ
    end

    @testset "0_10_100_50" begin
        dp = discretization_parameters(0, 10, 100)
        fp = frequency_parameters(dp, 50.)
        @test fp.ω == 50.
        @test fp.Ω == 250.
        @test fp.v ≈ v_0_10_100_50 atol = ϵ
        @test fp.K ≈ 1.2049415264262933 - 4.852640097709027im  atol = ϵ
    end

    @testset "0_5_100_10" begin
        dp = discretization_parameters(0, 5, 100)
        fp = frequency_parameters(dp, 10.)
        @test fp.ω == 10.
        @test fp.Ω == 25.
        @test fp.v ≈ v_0_5_100_10 atol = ϵ
        @test fp.K ≈ 2.478007029658684 - 0.33087937524443256im  atol = ϵ
    end

    @testset "0_10_200_10" begin
        dp = discretization_parameters(0, 10, 200)
        fp = frequency_parameters(dp, 10.)
        @test fp.ω == 10.
        @test fp.Ω == 50.
        @test fp.v ≈ v_0_10_200_10 atol = ϵ
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
        @test x == x_100_0_10
    end

    @testset "10, [0., 1., 10.], [100, 100]" begin
        x, fx, v = precompute_parameters(10, [0., 1., 10.], [100, 100])
        @test length(fx) == 202
        @test length(x) == 202
        @test length(v) == 202

        @test findmin(fx)[1] == 0 && findmax(fx)[1] == 0 
        @test x == x_100_100_0_1_10
    end
end

@testset "f_evolve_1!" begin
    @testset "1" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_1!(fx, x, Ct, 1)
        @test old_x == x
        fx == fx_evolve1_1
    end

    @testset "10" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_1!(fx, x, Ct, 10)
        @test old_x == x
        fx == fx_evolve1_10
    end
end

@testset "f_evolve_2!" begin
    @testset "1" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_2!(fx, x, 1)
        @test old_x == x
        fx == fx_evolve2_1
    end

    @testset "10" begin
        x, fx, v = precompute_parameters(1)
        old_x = copy(x)
        Ct = @. exp(-x^2 * 0.36)
        f_evolve_2!(fx, x, 10)
        @test old_x == x
        fx == fx_evolve2_10
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
