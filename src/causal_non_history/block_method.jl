
struct BlockMethod{T <: Number}
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
    sr_ζ::Vector{Vector{T}}
    sr_w::Vector{Vector{T}}
    sr_F::Vector{Vector{T}}
    sr_expt::Vector{Vector{T}}
    sr_expNout::Vector{Vector{T}}
    sr_Ic::Vector{T}
    sr_Icout::Vector{T}
    compute_self_response::Bool
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
    strength = 1.
end

function prepare_containers(setup::Setup, sources, ϵ, Nt, constants::Constants, containers=nothing; Q=1., compute_self_response=true)
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

    sr_ζ = Vector{Float64}[]
    sr_w = Vector{Float64}[]
    sr_F = Vector{Float64}[]
    sr_expt = Vector{Float64}[]
    sr_expNout = Vector{Float64}[]
    sr_Ic = Float64[]
    sr_Icout = Float64[]

    if compute_self_response
        for i in 1:Ns
            sr_setup = self_setup(setup, sources[i])
            precomp = precompute_parameters(sr_setup, params=constants)
            push!(sr_ζ, precomp.x)
            push!(sr_w, precomp.w)
            push!(sr_F, precomp.fx)
            push!(sr_expt, @. exp(-precomp.x^2*Δt̃))
            push!(sr_Ic, precomp.I_c)
            push!(sr_expNout,  @. exp(-precomp.x^2 * N[1] * Δt̃))
            push!(sr_Icout, constant_integral(sr_setup, constants, N[1]))
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
        qaux, 
        N,
        sr_ζ,
        sr_w,
        sr_F,
        sr_expt,
        sr_expNout,
        sr_Ic,
        sr_Icout,
        compute_self_response
    )
end

function evolve!(I, q, block::BlockMethod{T}) where {T <: Number}
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qaux, sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout, compute_self_response = block

    #if isempty(Kranges) return end

    Nb = size(q)[1]
    Nt = size(q)[2]
    K = size(load_delays)[1]
    for nt in 1:Nt
        current_q = map(x->x[1], block.load_buffer)
        @views for (q, b) in zip(q[:, nt], block.load_buffer)
            push!(b, q)
        end

        if compute_self_response
            for j in 1:Nb
                @. sr_F[j] = sr_expt[j] * (sr_F[j] - q[j, nt] / sr_ζ[j] + sr_expNout[j] * current_q[j] / sr_ζ[j])
                I[j, nt] += dot(sr_F[j], sr_w[j]) + q[j, nt] * sr_Ic[j] - current_q[j] * sr_Icout[j]
                @. sr_F[j] = sr_F[j] + (q[j, nt] - sr_expNout[j] * current_q[j]) / sr_ζ[j]
            end
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

        
        for target in 1:Nb
            for source in 1:Nb
                if !compute_self_response && source == target continue end
                @inbounds range = Kranges[K_min[source, target]]
                @inbounds @views I[target, nt] += dot(F[range, source], HM[range, target, source])
            end
        end
    end
end
