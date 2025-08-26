
struct BlockMethod{T <: Number}
    p::Vector{T}
    ζ::Vector{T}
    F::Matrix{T}
    expt::Vector{T}
    expNin::Vector{T}
    expNout::Vector{T}
    HM::Array{T, 3}
    load_delays::Matrix{CircularBuffer{T}}
    load_buffer::Vector{CircularBuffer{T}}
    ranges::Vector{UnitRange{Int}}
    Kranges::Vector{UnitRange{Int}}
    K_min::Matrix{Int}
    qaux::Vector{T} 
    N::Vector{Int}
    g::Vector{Vector{T}}
    compute_first_block::Bool
end

@with_kw struct PointSource{T <: Number} @deftype T
    x
    y
    z
    rb = 0.1
end
@with_kw struct LineSource{T <: Number} @deftype T
    x
    y
    D
    H
    rb = 0.1
end

function prepare_containers(setup::Setup, sources, ϵ, Nt, constants::Constants, containers=nothing; Q=1., compute_first_block=true)
    @unpack Δt, α, rb, kg, Δt̃ = constants

    n = 10
    expected_blocks = 5
    ϵ´ = 100ϵ/Q
    Ns = length(sources)
    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    distances = zeros(Ns, Ns)
    NR = zeros(Int, Ns, Ns)
    K_min = zeros(Int, Ns, Ns)

    for j in 1:Ns
        for i in 1:j-1
            distance, N_r = compute_distance(setup, sources[i], sources[j], constants, ϵ, Nt, Q)
            distances[i, j] = distance
            distances[j, i] = distance
            if N_r > Nt N_r = Nt end
            NR[i, j] = N_r
            NR[j, i] = N_r
        end
    end

    N, ND = choose_blocks(setup, sources, Nt, ϵ/expected_blocks, constants, Q)
    K = length(N) - 1

    for j in 1:Ns
        for i in 1:j-1
            last_block = findlast(x -> x <= NR[i, j], N)
            Km = min(isnothing(last_block) ? 1 : last_block, K)
            K_min[i, j] = Km
            K_min[j, i] = Km
        end
        K_min[j, j] = 1
    end
   
    ζ = zeros(0)
    W = zeros(0)
    indices = zeros(Int64, K+1)

    compute_ζ_discretization!(ζ, W, indices, setup; sources=sources, N=N, ND=ND, n=n, ϵ=ϵ/K, ϵ´=ϵ´, constants=constants, Q=Q)
    @views ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    @views Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    g = Vector{Float64}[]
    if compute_first_block 
        first_block_times = Δt .* (1:N[1])
        for i in 1:Ns
            sr_setup = self_setup(setup, sources[i])
            v = step_response.(first_block_times, Ref(sr_setup), Ref(constants); ϵ=ϵ´)
            if setup.image_strength != 0.
                sr_image_setup = image(sr_setup)
                v += setup.image_strength .* step_response.(first_block_times, Ref(sr_image_setup), Ref(constants); ϵ=ϵ´)
            end
            push!(g, v)
        end
    end

    # Preallocate objects
    F = zeros(length(ζ), Ns)
    expt = @. exp(-ζ^2*Δt̃)
    expNin = zeros(length(ζ))
    expNout = zeros(length(ζ))

    for i in 1:K
        @inbounds @. @views expNin[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i]*Δt̃)
        if i < K
            @inbounds @. @views expNout[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i+1]*Δt̃)
        end
    end

    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K, _ in 1:Ns]
    load_buffer = [CircularBuffer{Float64}(N[1]) for _ in 1:Ns]
    for buffer in load_buffer
        fill!(buffer, 0.)
    end
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    HM = zeros(length(ζ), length(sources), length(sources))

    compute_H!(HM, setup; ζ=ζ, W=W, expt=expt, sources=sources, distances=distances, constants=constants, ϵ=ϵ´, containers=containers, ND=ND, ranges=ranges)

    qaux = zeros(length(ζ))

    C = 1 / (2π^2*kg)

    BlockMethod(
        (1 .- expt) .* C .* W ./ζ, 
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
        qaux, 
        N,
        g,
        compute_first_block
    )
end

function evolve_F!(q,start,stop, block::BlockMethod{T}) where {T <: Number}
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qaux, g, compute_first_block = block

    #if isempty(Kranges) return end

    Nb = size(q)[1]
    Nt = size(q)[2]
    K = size(load_delays)[1]
    for nt in start:stop
        current_q = map(x->x[1], load_buffer)
        @views for (q, b) in zip(q[:, nt], load_buffer)
            push!(b, q)
        end

        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                push!(load_delays[i, j], qin)
                current_q[j] = qout
                @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
            end
            @views @. F[:, j] = expt * F[:, j] + qaux
        end

        if compute_first_block
            for i in 1:Nb
                Δq = diff([0; collect(load_buffer[i])])
                @inbounds I[i, nt] += dot(Δq, reverse(g[i]))
            end
        end
    end
end

function fmm_evaluation!(res,sources,targets,block::BlockMethod{T}) where {T <: Number}
        @unpack ζ, p, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qaux, g, compute_first_block = block

        # nζ,nt = length(ζ),length(targets)
    
        nζ = length(ζ)

        # targets = hcat([[p.x,p.y,p.z] for p in  positions]...)
        for k in 1:nζ
            zk = complex(block.ζ[k])
            charges = complex(block.F[k,:]) * p[k]
            vals = hfmm3d(1e-12,zk,sources,charges=charges,targets = targets, pg=1)
            @. res += imag(vals.pot)
        end
        
end



function evolve!(I, q, block::BlockMethod{T}) where {T <: Number}
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qaux, g, compute_first_block = block

    #if isempty(Kranges) return end

    Nb = size(q)[1]
    Nt = size(q)[2]
    K = size(load_delays)[1]
    for nt in 1:Nt
        current_q = map(x->x[1], load_buffer)
        @views for (q, b) in zip(q[:, nt], load_buffer)
            push!(b, q)
        end

        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                push!(load_delays[i, j], qin)
                current_q[j] = qout
                @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
            end
            @views @. F[:, j] = expt * F[:, j] + qaux
        end

        if compute_first_block
            for i in 1:Nb
                Δq = diff([0; collect(load_buffer[i])])
                @inbounds I[i, nt] += dot(Δq, reverse(g[i]))
            end
        end
  
        for target in 1:Nb
            for source in 1:Nb
                if !compute_first_block && source == target continue end
                if K_min[source, target] == 0 || length(Kranges) < K_min[source, target] continue end
                @inbounds range = Kranges[K_min[source, target]]
                @inbounds @views I[target, nt] += dot(F[range, source], HM[range, target, source])
            end
        end
    end
end

function build_nodes_and_weights!(ζ, W, segbuf, n)
    n_seg = length(segbuf)

    Nζ = n_seg*n
    Hs = Int(floor(n/2))
    p = n%2

    append!(ζ, zeros(Nζ))
    append!(W, zeros(Nζ))
    x, _, w = QuadGK.cachedrule(Float64, n)

    for (i, segment) in enumerate(segbuf)
        m = (segment.b-segment.a)/2
        c = (segment.b+segment.a)/2 
        @inbounds @views @. ζ[end-(n_seg-i+1)*n+1:end-(n_seg-i)*n-Hs] = m * x[2:2:end-1+p] + c
        @inbounds @views @. ζ[end-(n_seg-i+1)*n+1+Hs+p:end-(n_seg-i)*n] = -m * x[end-1-p:-2:2] + c
        @inbounds @views @. W[end-(n_seg-i+1)*n+1:end-(n_seg-i)*n-Hs] = m * w
        @inbounds @views @. W[end-(n_seg-i+1)*n+1+Hs+p:end-(n_seg-i)*n] = m * w[end-p:-1:1]
    end
end