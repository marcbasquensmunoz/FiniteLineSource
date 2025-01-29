using FiniteLineSource
using FiniteLineSource: create_bin_containers, get_bin
using FastGaussQuadrature
using LegendrePolynomials
using HMatrices
using DataStructures
using QuadGK
using SpecialFunctions
using Parameters

ϵ = 1e-6
Δt = 3600.
Nt = 900

bn = 10
bm = 10

α = 1e-6
kg = 3.
rb = 0.1
Δt̃ = Δt*α/rb^2

D = 0.
H = 100.
B = 10.
z_eval = D+H/2

q = [1. for t in 1:Nt]

bh_positions = [(B*(i-1)^2, B*(j-1)^2) for i in 1:bn for j in 1:bm]

compute_distance(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2)

params = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

struct BlockMethod{T <: Number}
    ζ::Vector{T}
    F::Vector{T}
    expt::Vector{T}
    expNin::Vector{T}
    expNout::Vector{T}
    HM::Array{T, 3}
    load_delays::Vector{Queue{T}}
    load_buffer::Queue{T}
    ranges::Vector{UnitRange{Int}}
    Kranges::Vector{UnitRange{Int}}
    K_min::Matrix{Int}
    qinaux::Vector{T}
    qoutaux::Vector{T}
end

function prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel)
    @unpack Δt, α, rb = params
    Δt̃ = Δt*α/rb^2

    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    R = filter!(e -> e != 0, [compute_distance(i, j) for i in bh_positions for j in bh_positions])
    unique!(R)

    # Choose blocks and compute points
    Nr = [compute_N_line(r, D, H, z_eval, ϵ, params) for r in R]
    sort!(Nr)
    N = choose_blocks(Nr, p = 10)

    ζ = [zeros(0) for _ in eachindex(N)]
    W = [zeros(0) for _ in eachindex(N)]

    for i in eachindex(N)
        Ni = N[i]
        No = i == length(N) ? 3*N[i] : N[i+1]
        ζ[i], W[i] = compute_ζ_points_line(Ni, No, ϵ, 1., 10, D, H, z_eval, params)
    end

    # Preallocate objects
    F = [zeros(length(ζ[i])) for i in eachindex(N)]
    expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(F)]
    expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(F)]
    expNout = [i == length(F) ? zeros(length(F[i])) : @. exp(-ζ[i]^2*(N[i+1])*Δt̃) for i in eachindex(F)]

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

    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    X, _, Fx, PP, XT = create_bin_containers(bins)

    distances = [compute_distance(source, target) for source in bh_positions, target in bh_positions]
    NR = [σ == 0 ? 0. : compute_N_line(σ, D, H, z_eval, ϵ, params) for σ in distances]
    K_min = [N_r == 0 ? 0 : findlast(x -> x <= N_r, N) for N_r in NR]
    indices = cumsum([length(F[i]) for i in eachindex(F)])
    insert!(indices, 1, 0)
    ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    ζflat = reduce(vcat, ζ)
    Wflat = reduce(vcat, W)
    Fflat = reduce(vcat, F)
    exptflat = reduce(vcat, expt)
    expNinflat = reduce(vcat, expNin)
    expNoutflat = reduce(vcat, expNout)

    HMflat = zeros(length(ζflat), length(bh_positions), length(bh_positions))

    C = 1 / (2π^2*kg)

    for (k, ζζ) in enumerate(ζflat), (j, target) in enumerate(bh_positions), i in 1:j-1
        σ = compute_distance(bh_positions[i], target)
        setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ, nmodel)
        bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
        bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
        HMflat[k, i, j] = C * Wflat[k] * compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=Fx[bin1], f2=Fx[bin2])
        HMflat[k, j, i] = HMflat[k, i, j]
    end

    qin = zeros(length(ζflat))
    qout = zeros(length(ζflat))

    N, BlockMethod(ζflat, Fflat, exptflat, expNinflat, expNoutflat, HMflat, load_delays, load_buffer, ranges, Kranges, K_min, qin, qout)
end

function evolve!(I, q, block::BlockMethod)
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, ranges, Kranges, K_min, qinaux, qoutaux = block

    for (nt, qt) in enumerate(q)
        enqueue!(load_buffer, qt)
        current_q = dequeue!(load_buffer)

        for i in eachindex(load_delays)
            qin = current_q
            enqueue!(load_delays[i], qin)
            qout = dequeue!(load_delays[i])
            current_q = qout
            @. qinaux[ranges[i]] = qin
            @. qoutaux[ranges[i]] = qout
        end
        @. F = expt * F + (qinaux * expNin - qoutaux * expNout) * (1 - expt) / ζ

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


N, block = @time prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel);
Istp = zeros(length(bh_positions), Nt)
@time evolve!(Istp, q, block)

#####################################
# Validation
#####################################

STP_test = zeros(length(bh_positions))
for target in eachindex(bh_positions)
    for source in eachindex(bh_positions)
        if source == target continue end
        σ = compute_distance(bh_positions[source], bh_positions[target])
        I, E = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z_eval)^2) / sqrt(4α * Δt * Nt))/(4*π*kg*sqrt(σ^2 + (zp - z_eval)^2)), D, D+H, atol = ϵ/100)
        STP_test[target] += I
    end
end

err = @. abs(Istp[:, end] - STP_test)

#####################################
# Validation of one case
#####################################
function compute_one_case(source, target)
    if source == target return 0. end
    σ = compute_distance(bh_positions[source], bh_positions[target])

    I_1 = 0.
    N_r = compute_N_line(σ, D, H, z_eval, ϵ, params)
    k_min = findlast(x -> x <= N_r, N)
    for j in length(F):-1:k_min
        for k in eachindex(ζ[j])
            I_1 += C * F[j][k] * W[j][k] * HM[j][target, source, k]
        end
    end

    I, E = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z_eval)^2) / sqrt(4α * Δt * Nt))/(4*π*kg*sqrt(σ^2 + (zp - z_eval)^2)), D, D+H, atol = ϵ/100)
    err = abs(I_1 - I)
end

cases = [compute_one_case(source, target) for target in 1:length(bh_positions), source in 1:length(bh_positions)]

compute_one_case(1, 4)
