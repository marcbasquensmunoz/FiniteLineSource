using FiniteLineSource
using Test

@testset "segment_to_segment_step_response" begin
    source_segment = BoreholeSegment(0, 0, 1, 5, 0.5)
    receptor_segment = BoreholeSegment(0, 1, 1, 5, 0.5)
    @test segment_to_segment_step_response(0, source=source_segment, receptor=receptor_segment) == 0
    @test segment_to_segment_step_response(1000, source=source_segment, receptor=receptor_segment) < 0
    @test segment_to_segment_step_response(1000, source=source_segment, receptor=receptor_segment)     ≈ 0        atol = 5*10^-10
    @test segment_to_segment_step_response(100000, source=source_segment, receptor=receptor_segment)   ≈ -0.04414 atol = 5*10^-6
    @test segment_to_segment_step_response(10000000, source=source_segment, receptor=receptor_segment) ≈ -0.12464 atol = 5*10^-6
end
