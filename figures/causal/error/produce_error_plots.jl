using FiniteLineSource

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q_step(t) = 1.
q_synth(t) = 20*sin(2π*t/8760) + 5*sin(2π*t/24) + 5.

r̃_range = 10 .^ collect(1:1/6:3)
ϵ_range = 10. .^ collect(-2:-2:-12)

include("p2p_error.jl")
include("l2p_error.jl")
include("l2l_error.jl")