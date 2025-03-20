
struct BlockMethod{T <: Number}
    ζ::Vector{T}
    F::Vector{T}
    expt::Vector{T}
    expNin::Vector{T}
    expNout::Vector{T}
    HM::Array{T, 3}
    load_delays::Vector{CircularBuffer{T}}
    load_buffer::CircularBuffer{T}
    ranges::Vector{UnitRange{Int}}
    Kranges::Vector{UnitRange{Int}}
    K_min::Matrix{Int}
    qinaux::Vector{T}
    qoutaux::Vector{T}
    N::Vector{Int}
end
@with_kw struct LineSource{T <: Number} @deftype T
    x
    y
    D
    H
end

function prepare_containers(setup::Setup, sources, ϵ, Nt, constants::Constants, containers=nothing)
    @unpack Δt, α, rb, kg, Δt̃ = constants

    n = 10
    ϵ´ = 100ϵ
    Ns = length(sources)
    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    distances = zeros(Ns, Ns)
    NR = zeros(Int, Ns, Ns)
    K_min = zeros(Int, Ns, Ns)

    for j in 1:Ns
        for i in 1:j-1
            distance, N_r = compute_distance(setup, sources[i], sources[j], constants, ϵ, Nt)
            distances[i, j] = distance
            distances[j, i] = distance
            if N_r > Nt N_r = Nt end
            NR[i, j] = N_r
            NR[j, i] = N_r
        end
    end

    Nr = filter!(e -> e != 0, unique(NR))
    N, ND = choose_blocks(setup, sources, distances, Nr, Nt, ϵ/5, constants)
    K = length(N) - 1

    #@show N, ND
    for j in 1:Ns
        for i in 1:j-1
            last_block = findlast(x -> x <= NR[i, j], N)
            Km = min(isnothing(last_block) ? 1 : last_block, K)
            K_min[i, j] = Km
            K_min[j, i] = Km
        end
    end
   
    ζ = zeros(0)
    W = zeros(0)
    indices = zeros(Int64, K+1)

    compute_ζ_discretization!(ζ, W, indices, setup; sources=sources, N=N, ND=ND, n=n, ϵ=ϵ/K, ϵ´=ϵ´, constants=constants, containers=containers)
    @views ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    @views Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    #sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout = bakhalov_discretization(10, setup)
    
    # Preallocate objects
    F = zeros(length(ζ))
    expt = @. exp(-ζ^2*Δt̃)
    expNin = zeros(length(ζ))
    expNout = zeros(length(ζ))

    for i in 1:K
        @inbounds @. @views expNin[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i]*Δt̃)
        if i < K
            @inbounds @. @views expNout[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i+1]*Δt̃)
        end
    end
    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K]
    load_buffer = CircularBuffer{Float64}(N[1])
    fill!(load_buffer, 0.)
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    HM = zeros(length(ζ), length(sources), length(sources))

    compute_H!(HM, setup; ζ=ζ, W=W, expt=expt, sources=sources, distances=distances, constants=constants, ϵ=ϵ´, containers=containers, N=N, ND=ND, ranges=ranges)

    qin = zeros(length(ζ))
    qout = zeros(length(ζ))

    BlockMethod(
        ζ, 
        F, 
        expt, 
        expNin, 
        expNout, 
        HM, 
        load_delays, 
        load_buffer, 
        ranges, 
        Kranges, 
        K_min, 
        qin, 
        qout,
        N
        #=sr_ζ,
        sr_w,
        sr_F,
        sr_expt,
        sr_expNout,
        sr_Ic,
        sr_Icout=#
    )
end

function evolve!(I, q, block::BlockMethod{T}) where {T <: Number}
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qinaux, qoutaux#=, sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout =#= block

    if isempty(Kranges) return end

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
