using FiniteLineSource
using Test
using SpecialFunctions

include("Aqua.jl")

const ϵ = 5*10^-14
@testset "segment_to_segment_step_response" begin
    source_segment = BoreholeSegment(0, 0, 1, 5, 0.5)
    receptor_segment = BoreholeSegment(0, 1, 1, 5, 0.5)
    @test segment_to_segment_step_response(0, source=source_segment, receptor=receptor_segment)        ≈ 0        atol = 5*10^-20
    @test segment_to_segment_step_response(100, source=source_segment, receptor=receptor_segment)      ≈ 0        atol = 5*10^-20
    @test segment_to_segment_step_response(1000, source=source_segment, receptor=receptor_segment)     ≈ 0        atol = 5*10^-10
    @test segment_to_segment_step_response(100000, source=source_segment, receptor=receptor_segment)   ≈ 0.131890 atol = 5*10^-6
    @test segment_to_segment_step_response(10000000, source=source_segment, receptor=receptor_segment) ≈ 1.074141 atol = 5*10^-6
end

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