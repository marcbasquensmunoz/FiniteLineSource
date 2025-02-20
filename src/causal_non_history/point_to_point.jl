
compute_distance_3D(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)

"""
Compute the number of steps N that can be skipped for a given distance r
"""
function compute_N(r, ϵ, params::Constants) 
    @unpack Δt, α, kg = params
    Int(floor((r / sqrt(4α) / erfcinv(4*π*r*kg*ϵ))^2/Δt))
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points!(ζ, W, N, No, ϵ, Q, n, params::Constants)
    @unpack Δt, α, rb = params
    Δt̃ = Δt*α/rb^2

    heatwave(r) = erfc(r/sqrt(4α*N*Δt)) / r - ϵ
    r = find_zero(heatwave, sqrt(N))
    a = 0.
    b = sqrt(-log(ϵ/Q) / (N*Δt̃))
    r̃ = r/rb

    guide(ζ) = (exp(-ζ^2*N*Δt̃) + 100 * exp(-ζ^2*No*Δt̃)) * sin(r̃*ζ) / (r*ζ) * (1 - exp(-ζ^2*Δt̃))
    #guide(ζ) = exp(-ζ^2*N*Δt̃) * sin(r̃*ζ) / (r*ζ) * (1 -  exp(-ζ^2*Δt̃))
    _, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    sort!(segbuf, by=s->s.a)
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
        @views @. ζ[end-(n_seg-i+1)*n+1:end-(n_seg-i)*n-Hs] = m * x[2:2:end-1+p] + c
        @views @. ζ[end-(n_seg-i+1)*n+1+Hs+p:end-(n_seg-i)*n] = -m * x[end-1-p:-2:2] + c
        @views @. W[end-(n_seg-i+1)*n+1:end-(n_seg-i)*n-Hs] = m * w
        @views @. W[end-(n_seg-i+1)*n+1+Hs+p:end-(n_seg-i)*n] = m * w[end-p:-1:1]
    end
    return Nζ
end

function prepare_containers_ptp(positions, ϵ, Nt, params)
    @unpack Δt, α, rb, kg = params
    Δt̃ = Δt*α/rb^2
    n = 10

    Ns = length(positions)
    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    distances = zeros(Ns, Ns)
    NR = zeros(Int, Ns, Ns)
    K_min = zeros(Int, Ns, Ns)

    for j in 1:Ns
        for i in 1:j-1
            r = compute_distance_3D(positions[i], positions[j])
            N_r = compute_N(r, ϵ, params)
            distances[i, j] = r
            distances[j, i] = r
            NR[i, j] = N_r
            NR[j, i] = N_r
        end
    end

    Nr = filter!(e -> e != 0, unique(NR))
    N = choose_blocks(Nr, Nt, p = 10)
    K = length(N) - 1

    for j in 1:Ns
        for i in 1:j-1
            Km = findlast(x -> x <= NR[i, j], N)
            K_min[i, j] = Km
            K_min[j, i] = Km
        end
    end

    ζ = zeros(0)
    W = zeros(0)
    indices = zeros(Int64, K+1)

    for i in 1:K
        N_block = compute_ζ_points!(ζ, W, N[i], N[i+1], ϵ, 1., n, params)
        @views indices[i+1:end] .+= N_block
    end

    @views ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    @views Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    # Preallocate objects
    F = zeros(length(ζ))
    expt = @. exp(-ζ^2*Δt̃)
    expNin = zeros(length(ζ))
    expNout = zeros(length(ζ))
    for i in 1:K
        @. @views expNin[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i]*Δt̃)
        if i < K
            @. @views expNout[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i+1]*Δt̃)
        end
    end

    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K]
    load_buffer = CircularBuffer{Float64}(N[1])
    fill!(load_buffer, 0.)
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    HM = zeros(length(ζ), length(positions), length(positions))

    C = 1 / (2π^2*kg)

    for (k, ζζ) in enumerate(ζ), j in eachindex(positions), i in 1:j-1
        r = distances[i, j]
        HM[k, i, j] = C * W[k] * sin(r/rb*ζζ) / (r*ζζ) * (1 - expt[k])
        HM[k, j, i] = HM[k, i, j]
    end

    qin = zeros(length(ζ))
    qout = zeros(length(ζ))

    N, BlockMethod(ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, ranges, Kranges, K_min, qin, qout)
end
    
function evolve_ptp!(I, q, block::BlockMethod)
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, ranges, Kranges, K_min, qinaux, qoutaux = block

    for (nt, qt) in enumerate(q)
        current_q = load_buffer[1]
        push!(load_buffer, qt)

        for i in eachindex(load_delays)
            qin = current_q
            qout = load_delays[i][1]
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