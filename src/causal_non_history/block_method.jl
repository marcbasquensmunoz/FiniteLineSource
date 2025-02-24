
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
end

"""
Computes the blocks to be used
"""
function choose_blocks(Nr, Nt; p = 10)
    Nmin = minimum(Nr)
    Ncurrent = Nmin
    N = zeros(Int, 0)

    while Ncurrent < Nt
        push!(N, Int(floor(Ncurrent)))
        Ncurrent *= p
    end
    !(Nmin in N) && push!(N, Nmin)
    sort!(N)
    push!(N, Nt)

    @views for (na, nb) in zip(N[1:end-1], N[2:end])
        if isempty(filter(n -> n in na:nb, Nr)) 
            deleteat!(N, findfirst(x->x==na, N))
        end
    end
    return N
end

"""
Compute the minimum number of terms n to use in the Legendre expansion of the function 1/sqrt(x^2-σ^2)
to obtain an error lower than ϵ in the interval [a, b].
"""
function N_bound(ϵ, a, b, σ, nmodel)
    Int(ceil(eval(nmodel, ϵ, σ, a, b)))
end

@with_kw struct LineKernelContainers{T <: Number} @deftype T
    x::Vector{T}
    X::Vector{T}
    f::Vector{T}
    P::Matrix{T}
end 

function create_bin_containers(bins)
    X = Vector{Vector{Float64}}()
    PP = Vector{Matrix{Float64}}()
    F = Vector{Vector{Float64}}()
    XT = Vector{Vector{Float64}}()

    for n in bins
        x, w = gausslegendre(n+1)
        P = zeros(n+1, n+1)
        for s in 1:n+1
            @inbounds @views collectPl!(P[:, s], x[s], lmax=n)
            @inbounds @views @. P[:, s] *= w[s]
        end
        f = zeros(n+1)
        xt = zeros(n+1)

        push!(X, x)
        push!(F, f)
        push!(PP, P)
        push!(XT, xt)
    end

    X, F, PP, XT
end
get_bin(n, bins) = findfirst(m -> n < m, bins)

function evolve!(I, q, block::BlockMethod{T}) where {T <: Number}
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
