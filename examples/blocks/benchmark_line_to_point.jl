using FiniteLineSource
using FiniteLineSource: create_bin_containers, get_bin, DiscretizationParameters, SegmentToPoint, BlockMethod
using FastGaussQuadrature
using LegendrePolynomials
using HMatrices
using DataStructures
using QuadGK
using SpecialFunctions
using Parameters
using LinearAlgebra
using Roots

ϵ = 1e-6
Δt = 3600.
Nt = 10000

bn = 1
bm = 2

α = 1e-6
kg = 3.
rb = 0.1
Δt̃ = Δt*α/rb^2

D = 0.
H = 100.
B = 1.1
z_eval = D+H/2

q = [1. for t in 1:Nt]

bh_positions = [(B*(i-1)^2, B*(j-1)^2) for i in 1:bn for j in 1:bm]

compute_distance(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2)

params = Constants(Δt=Δt, α=α, kg=kg, rb=rb,line_points=[1, 1, 1, 1, 1] .* 500, line_limits=[0., 0.1, 0.3, 0.7, 0.9, 1.])

function bakhalov_discretization(N, setup)
    @unpack D, H, z = setup
    σ = rb
    Δt̃ = α*Δt/rb^2
    n = 10

    z_int = log((z-D + sqrt(σ^2 + (z-D)^2)) / (z-D-H + sqrt(σ^2 + (z-D-H)^2)))
    a = 0.
    b = find_zero(bb -> 2ϵ/(z_int*rb) - N * (gamma(0, bb^2*Δt̃) - gamma(0, bb^2*2*Δt̃)), -log(ϵ) / (2N*Δt̃))

    guide(ζ) = (1 + exp(-ζ^2*Δt̃*N))* (1 - exp(-ζ^2*Δt̃)) / ζ
    _, _, segments = quadgk_segbuf(guide, a, b)
    dps = @views [DiscretizationParameters(s.a, s.b, n) for s in segments]
    x  = reduce(vcat, (dp.x for dp in dps))
    w  = reduce(vcat, [FiniteLineSource.precompute_coefficients(setup, dp=dp, params=params, containers=FiniteLineSource.EmptyContainer()) for (i, dp) in enumerate(dps)])
    fx = zeros(sum([dp.n+1 for dp in dps]))
    perm = sortperm(x)

    expt = @. exp(-x^2 * Δt̃)
    exptout = @. exp(-x^2 * N * Δt̃)
    Ic = log((z-D + sqrt(σ^2 + (z-D)^2))/(z-D-H + sqrt(σ^2 + (z-D-H)^2))) /  (4π * kg)
    Icout = quadgk(zp -> erf(sqrt(rb^2 + (zp - z)^2)/rb/sqrt(4*N*Δt̃)) / sqrt(rb^2 + (zp - z)^2), D, D+H)[1] / (4π * kg)

    x[perm], w[perm], fx, expt, exptout, Ic, Icout
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
    N = isempty(Nr) ? [10] :  choose_blocks(Nr, Nt, p = 10)
    unique!(N)
    K = length(N) - 1

    ζ = [zeros(0) for _ in eachindex(N)]
    W = [zeros(0) for _ in eachindex(N)]

    for i in eachindex(N)
        Ni = N[i]
        No = i == length(N) ? 3*N[i] : N[i+1]
        ζ[i], W[i] = compute_ζ_points_line(Ni, No, ϵ, 1., 10, D, H, z_eval, params)
    end

    setup = SegmentToPoint(σ=rb, D=D, H=H, z=z_eval)
    #sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout = bakhalov_discretization(10, setup)
    
    # Preallocate objects
    F = [zeros(length(ζ[i])) for i in eachindex(N)]
    expt = [@. exp(-ζ[i]^2*Δt̃) for i in eachindex(F)]
    expNin = [@. exp(-ζ[i]^2*N[i]*Δt̃) for i in eachindex(F)]
    expNout = [i == length(F) ? zeros(length(F[i])) : @. exp(-ζ[i]^2*(N[i+1])*Δt̃) for i in eachindex(F)]


    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K]
    load_buffer = CircularBuffer{Float64}(N[1])
    fill!(load_buffer, 0.)
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    X, _, Fx, PP, XT = create_bin_containers(bins)

    distances = [compute_distance(source, target) for source in bh_positions, target in bh_positions]
    NR = [σ == 0 ? 0. : compute_N_line(σ, D, H, z_eval, ϵ, params) for σ in distances]
    K_min = [N_r == 0 ? 1 : findlast(x -> x <= N_r, N) for N_r in NR]
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

    for (k, ζζ) in enumerate(ζflat), (j, target) in enumerate(bh_positions), i in 1:j
        σ = i == j ? rb : compute_distance(bh_positions[i], target)
        setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ, nmodel)
        bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
        bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
        HMflat[k, i, j] = C * Wflat[k] * (1-exptflat[k]) / ζflat[k] * compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=Fx[bin1], f2=Fx[bin2])
        HMflat[k, j, i] = HMflat[k, i, j]
    end

    qin = zeros(length(ζflat))
    qout = zeros(length(ζflat))

    BlockMethod(
        ζflat, 
        Fflat, 
        exptflat, 
        expNinflat, 
        expNoutflat, 
        HMflat, 
        load_delays, 
        load_buffer, 
        ranges, 
        Kranges, 
        K_min, 
        qin, 
        qout,
        #=sr_ζ,
        sr_w,
        sr_F,
        sr_expt,
        sr_expNout,
        sr_Ic,
        sr_Icout=#
    )
end


function evolve!(I, q, block::BlockMethod)
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qinaux, qoutaux#=, sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout =#= block

    for (nt, qt) in enumerate(q)
        current_q = load_buffer[1]
        push!(load_buffer, qt)
        
        #=
        @. sr_F = sr_expt * (sr_F - qt / sr_ζ + sr_expNout * current_q / sr_ζ)
        for target in 1:size(block.K_min)[1]
            I[target, nt] += dot(sr_F, sr_w) + qt * sr_Ic - current_q * sr_Icout
        end
        @. sr_F = sr_F + (qt - sr_expNout * current_q) / sr_ζ
        =#
        for i in eachindex(load_delays)
            qin = current_q
            qout = i == length(load_delays) ? 0. : load_delays[i][1]
            push!(load_delays[i], qin)
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

nmodel = FiniteLineSource.load_nmodel()
block = @time prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel);
Istp = zeros(length(bh_positions), Nt)
@time evolve!(Istp, q, block)

#####################################
# Validation
#####################################

STP_test = zeros(length(bh_positions))
for target in eachindex(bh_positions)
    for source in eachindex(bh_positions)
        σ = source == target ? rb : compute_distance(bh_positions[source], bh_positions[target])
        if source == target continue end
        I, E = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z_eval)^2) / sqrt(4α * Δt * Nt))/(4*π*kg*sqrt(σ^2 + (zp - z_eval)^2)), D, D+H, atol = ϵ/100)
        STP_test[target] += I
    end
end

err = @. abs(Istp[:, end] - STP_test)

#=
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
=#

Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=0.01 prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel);
PProf.Allocs.pprof(from_c=false)