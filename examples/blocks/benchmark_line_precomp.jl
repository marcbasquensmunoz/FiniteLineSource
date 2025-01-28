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
Nt = 24*8760

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

    for j in 1:N[1]
        enqueue!(load_buffer, 0.)
    end
    for i in 1:length(load_delays)-1
        for j in 1:N[i+1]-N[i]
            enqueue!(load_delays[i], 0.)
        end
    end

    #HM = Vector{Vector{HMatrix{ClusterTree{2, Float64}, Float64}}}[]
    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    X, _, Fx, PP, XT = create_bin_containers(bins)

    #=
    function G(x, y, k)
        if x == y return 1. end
        σ = sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
        setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=σ)
        lineparams = LineKernelParams(setup, k/rb, ϵ, eval)
        bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
        bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
        compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=Fx[bin1], f2=Fx[bin2])
    end=#
    #allζ = reduce(vcat, ζ)
    #HM = [zeros(length(bh_positions), length(bh_positions), length(ζζ)) for ζζ in ζ]
    HM = [[zeros(length(bh_positions), length(bh_positions)) for ζζ in ζ[i]] for i in eachindex(ζ)]

    for (l, vζ) in enumerate(ζ), (k, ζζ) in enumerate(vζ), (j, target) in enumerate(bh_positions), i in 1:j-1
        σ = compute_distance(bh_positions[i], target)
        setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ, nmodel)
        bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
        bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
        HM[l][k][i, j] = compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=Fx[bin1], f2=Fx[bin2])
        HM[l][k][j, i] = HM[l][k][i, j]
    end

    N, ζ, W, F, expt, expNin, expNout, HM, load_delays, load_buffer
end

N, ζ, W, F, expt, expNin, expNout, HM, load_delays, load_buffer = @time prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel)
@info "Precomputation finished"
#=
Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=0.05 prepare_containers(bh_positions, D, H, z_eval, ϵ, params, nmodel)
PProf.Allocs.pprof(from_c=false)
=#

@time begin
#####################################
# Evolution in time
#####################################
for qt in q
    enqueue!(load_buffer, qt)
    current_q = dequeue!(load_buffer)

    for i in eachindex(F)
        qin = current_q
        enqueue!(load_delays[i], qin)
        qout = dequeue!(load_delays[i])
        current_q = qout
        @. F[i] = expt[i] * F[i] + qin * expNin[i] - qout * expNout[i]
    end
end

#####################################
# Complete F
#####################################
C = 1 / (2π^2*kg)
for i in eachindex(F)
    @. F[i] = F[i] * (1 - expt[i]) / ζ[i]
end

#####################################
# Compute the integral
#####################################

Istp = zeros(length(bh_positions))
for target in eachindex(bh_positions)
    for source in eachindex(bh_positions)
        if source == target continue end
        r = compute_distance(bh_positions[source], bh_positions[target])
        N_r = compute_N_line(r, D, H, z_eval, ϵ, params)
        k_min = findlast(x -> x <= N_r, N)
        for j in length(F):-1:k_min
            for k in eachindex(ζ[j])
                Istp[target] += C * F[j][k] * W[j][k] * HM[j][k][target, source]
            end
        end
    end
end
end

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

err = @. abs(Istp - STP_test)

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

##################
# Benchmarking
##################
# Here I use H-matrices for the computation and storage of the kernels.
# It is not so useful in terms of total memory size, but it appears to be faster
@time begin
for i in eachindex(ζ)
    for ζζ in ζ[i]

        K = KernelMatrix{Function, Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}, Float64}(bh_positions, bh_positions) do x, y
            if x == y return 1. end
            setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2))
            lineparams = LineKernelParams(setup, ζζ/rb, ϵ, eval)

            bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
            bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)

            compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=F[bin1], f2=F[bin2])
        end
    
        assemble_hmatrix(K, atol=ϵ)
    end
end
end

# Here I store the matrices myself
@time begin
function G(x, y, k)
    if x == y return 1. end
    setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2))
    lineparams = LineKernelParams(setup, k/rb, ϵ, eval)
    bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
    bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
    compute_kernel_line(lineparams; x1=X[bin1], x2=X[bin2], X1=XT[bin1], X2=XT[bin2], P1=PP[bin1], P2=PP[bin2], f1=F[bin1], f2=F[bin2])
end
HM = [[[G(source, target, ζζ) for source in bh_positions, target in bh_positions] for ζζ in ζ[i]] for i in eachindex(ζ)]
nothing
end
