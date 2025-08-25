using FiniteLineSource
using FiniteLineSource: LineSource
using BenchmarkTools


ϵ = 1e-6
Nt = 8760*20
B = 1.

D1 = 0.
D2 = 0.
H1 = 100.
H2 = 100.

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)


q = hcat(
    [sin(t/8760) for t in 1:Nt],
    [sin(t/8760) for t in 1:Nt]
)'

bh_positions = [LineSource(x=0., y=0., D=D1, H=H1, rb=rb), LineSource(x=B, y=0., D=D2, H=H2, rb=rb)]

#####################################
# Error analysis
#####################################

setup = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=B, image_strength = -1.)
image_setup = SegmentToSegment(D1=-D1-H1, H1=H1, D2=D2, H2=H2, σ=B)
sr_setup = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=rb)
sr_image_setup = SegmentToSegment(D1=-D1-H1, H1=H1, D2=D2, H2=H2, σ=rb)

# Block method
containers = FiniteLineSource.AsymptoticContainers(10)
block = @time prepare_containers(setup, bh_positions, ϵ/length(bh_positions), Nt, constants, containers, Q=maximum(q), compute_first_block=true);
Ib = zeros(length(bh_positions), Nt)
@time evolve!(Ib, q, block)

# Convolution
C_2to1 = @time convolve_step(q[2,:], SegmentToSegmentOld(setup); params=constants, ϵ=ϵ/100)
C_2to1_image = @time convolve_step(-q[2,:], SegmentToSegmentOld(image_setup); params=constants, ϵ=ϵ/100)
C_1sr = @time convolve_step(q[1,:], SegmentToSegmentOld(sr_setup); params=constants, ϵ=ϵ/100)
C_1sr_image = @time convolve_step(-q[1,:], SegmentToSegmentOld(sr_image_setup); params=constants, ϵ=ϵ/100)

C1 = C_2to1 + C_2to1_image + C_1sr + C_1sr_image

# Error
err = @. abs(Ib[1, :] - C1)
maximum(err)
