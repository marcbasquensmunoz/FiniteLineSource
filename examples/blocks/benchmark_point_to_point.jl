using FiniteLineSource: create_bin_containers, get_bin, Constants, compute_N, choose_blocks, compute_ζ_points, BlockMethod,precompute_parameters, compute_integral_throught_history!, PointToPoint
using FastGaussQuadrature
using LegendrePolynomials
using HMatrices
using DataStructures
using QuadGK
using SpecialFunctions
using Parameters
using LinearAlgebra

ϵ = 1e-6
Δt = 3600.
Nt = 1000

bn = 2
bm = 1
bl = 1

α = 1e-6
kg = 3.
rb = 0.1
Δt̃ = Δt*α/rb^2

B = 5.

q = [1. for t in 1:Nt]

positions = [(B*(i-1)^2, B*(j-1)^2, B*(k-1)^2) for i in 1:bn for j in 1:bm for k in 1:bl]

compute_distance_3D(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)
params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

function prepare_containers_ptp(positions, ϵ, Nt, params)
    @unpack Δt, α, rb, kg = params
    Δt̃ = Δt*α/rb^2
    n = 10

    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    distances = [compute_distance_3D(source, target) for source in positions, target in positions]
    R = filter!(e -> e != 0, [compute_distance_3D(i, j) for i in positions for j in positions])
    unique!(R)

    # Choose blocks and compute points
    Nr = [compute_N(r, ϵ, params) for r in R]
    sort!(Nr)
    N = choose_blocks(Nr, Nt, p = 10)
    K = length(N) - 1

    ζ = [zeros(0) for _ in 1:K]
    W = [zeros(0) for _ in 1:K]

    x, w = gausslegendre(n)
    for i in 1:K
        ζ[i], W[i] = compute_ζ_points(N[i], N[i+1], ϵ, 1., n, params, x, w)
    end
    # Preallocate objects
    F = [zeros(length(ζ[i])) for i in 1:K]
    expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(F)]
    expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(F)]
    expNout = [i == length(F) ? zeros(length(F[i])) : @. exp(-ζ[i]^2*(N[i+1])*Δt̃) for i in eachindex(F)]

    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K]
    load_buffer = CircularBuffer{Float64}(N[1])
    fill!(load_buffer, 0.)
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    #=
    load_delays = [Queue{Float64}() for _ in eachindex(F)]
    load_buffer = Queue{Float64}()

    for _ in 1:N[1]
        enqueue!(load_buffer, 0.)
    end
    for i in 1:length(load_delays)-1
        for _ in 1:N[i+1]-N[i]
            enqueue!(load_delays[i], 0.)
        end
    end
    =#

    NR = [r == 0 ? 0 : compute_N(r, ϵ, params) for r in distances]
    K_min = [N_r == 0 ? 0 : findlast(x -> x <= N_r, N) for N_r in NR]
    indices = cumsum([length(F[i]) for i in eachindex(F)])
    insert!(indices, 1, 0)
    @views ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    @views Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    ζflat = reduce(vcat, ζ)
    Wflat = reduce(vcat, W)
    Fflat = reduce(vcat, F)
    exptflat = reduce(vcat, expt)
    expNinflat = reduce(vcat, expNin)
    expNoutflat = reduce(vcat, expNout)

    HMflat = zeros(length(ζflat), length(positions), length(positions))

    C = 1 / (2π^2*kg)

    for (k, ζζ) in enumerate(ζflat), (j, target) in enumerate(positions), i in 1:j-1
        r = compute_distance_3D(positions[i], target)
        HMflat[k, i, j] = C * Wflat[k] * sin(r/rb*ζζ) / (r*ζζ) * (1 - exptflat[k])
        HMflat[k, j, i] = HMflat[k, i, j]
    end

    qin = zeros(length(ζflat))
    qout = zeros(length(ζflat))

    N, BlockMethod(ζflat, Fflat, exptflat, expNinflat, expNoutflat, HMflat, load_delays, load_buffer, ranges, Kranges, K_min, qin, qout)
end

function evolve!(I, q, block::BlockMethod)
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, ranges, Kranges, K_min, qinaux, qoutaux = block

    for (nt, qt) in enumerate(q)
        #enqueue!(load_buffer, qt)
        #current_q = dequeue!(load_buffer)

        current_q = load_buffer[1]
        push!(load_buffer, qt)

        for i in eachindex(load_delays)
            qin = current_q
            qout = load_delays[i][1]
            push!(load_delays[i], qin)
            #enqueue!(load_delays[i], qin)
            #qout = dequeue!(load_delays[i])
            current_q = qout
            @. qinaux[ranges[i]] = qin
            @. qoutaux[ranges[i]] = qout
        end
        @. F = expt * F + (qinaux * expNin - qoutaux * expNout)

        bh_indices = 1:size(block.K_min)[1]
        for target in bh_indices
            for source in bh_indices
                if source == target continue end
                range = Kranges[K_min[source, target]]
                @views I[target, nt] += dot(F[range], HM[range, target, source])
            end
        end
    end
end


N, block = @time prepare_containers_ptp(positions, ϵ, Nt, params);
Iptp = zeros(length(positions), Nt)
@time evolve!(Iptp, q, block)

#####################################
# Validation
#####################################

#=
for target in eachindex(positions)
    for source in eachindex(positions)
        if source == target continue end
        r = compute_distance_3D(positions[source], positions[target])
        test[target] += erfc(r / sqrt(4α * Δt * Nt))/(4*π*kg*r)
    end
end
err = @. abs(Iptp[:, end] - test)
=#

I = zeros(Nt)
setup = PointToPoint(r = B)
precomp = @time precompute_parameters(setup, params=params)
@time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

err = @. abs(Iptp[1, :] - I)


Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=1  prepare_containers_ptp(positions, ϵ, Nt, params);
PProf.Allocs.pprof(from_c=false)

