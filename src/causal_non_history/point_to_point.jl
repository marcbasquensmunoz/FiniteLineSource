
compute_distance_3D(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)

"""
Compute the number of steps N that can be skipped for a given distance r
"""
function compute_N(r, ϵ, params::Constants, Nt=Inf) 
    @unpack Δt, α, kg = params
    if 4*π*r*kg*ϵ > 1
        return Nt
    end
    Int(floor((r / sqrt(4α) / erfcinv(4*π*r*kg*ϵ))^2/Δt))
end
function compute_r(N, ϵ, params::Constants) 
    @unpack Δt, α, kg = params
    f(r) = erfc(r / sqrt(4α*Δt*N)) / (4*π*r*kg) - ϵ
    problem = ZeroProblem(f, (0, N))
    solve(problem, Roots.Brent(), xatol = 1e-2)
end

function get_representative_ptp(sources)
    min_dist = Inf
    setup = PointToPoint(r=0.)
    for (i, s) in enumerate(sources)
        for t in @views sources[1:i-1]
            dist = compute_distance_3D(s, t)
            if dist < min_dist
                min_dist = dist
                setup = PointToPoint(r=dist)
            end
        end
    end
    return min_dist
end

"""
Computes the blocks to be used
"""
function choose_blocks(::PointToPoint, sources, Nt, ϵ, constants)
    Ncurrent = compute_N(get_representative_ptp(sources), ϵ, constants)
    N = Int[]
    while Ncurrent < Nt
        push!(N, Ncurrent)
        Ncurrent *= 10
    end
    push!(N, Nt)
    return N, map(n -> compute_r(n, ϵ, constants), N)
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points!(ζ, W, N, No, ϵ, n, params::Constants)
    @unpack Δt, α, rb, Δt̃, kg = params

    r = compute_r(N, ϵ, params) 
    a = 0.
    b = sqrt(-log(ϵ) / (N*Δt̃))
    r̃ = r/rb

    guide(ζ) = (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * sin(r̃*ζ) / (r*ζ) * (1 - exp(-ζ^2*Δt̃))
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

function compute_distance(::PointToPoint, source, target, params, ϵ, Nt)
    r = compute_distance_3D(source, target)
    N_r = compute_N(r, ϵ, params, Nt)
    r, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::PointToPoint; sources, N, ND, n, ϵ, ϵ´, constants)
    K = length(N) - 1
    for i in 1:K
        N_block = compute_ζ_points!(ζ, W, N[i], N[i+1], ϵ, n, constants)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::PointToPoint; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers, ND, ranges)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    for (k, ζζ) in enumerate(ζ), j in eachindex(sources), i in 1:j-1
        r = distances[i, j]
        HM[k, i, j] = C * W[k] * sin(r/rb*ζζ) / (r*ζζ) * (1 - expt[k])
        HM[k, j, i] = HM[k, i, j]
    end
end
