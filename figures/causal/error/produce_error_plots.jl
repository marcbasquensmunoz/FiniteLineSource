using FiniteLineSource

α = 1e-6
kg = 3.
rb = 0.1
Δt = 3600.
Nt = 8760 * 4

constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

q_step = ones(Nt)
q_synth = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:Nt]

r_range = 1.:1:50.
ϵ_range = 10. .^ collect(-2:-2:-12)

include("p2p_error.jl")
include("l2p_error.jl")
include("l2l_error.jl")