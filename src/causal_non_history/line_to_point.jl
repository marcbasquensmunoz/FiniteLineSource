using Plots 
compute_distance_2D(s, t) = sqrt((s.x - t.x)^2 + (s.y - t.y)^2)

@with_kw struct LineKernelParams{T <: Number} @deftype T
    r1
    r2
    r3
    ω
    σ
    ϵ
end

function LineKernelParams(setup::SegmentToPoint, ω, ϵ)
    @unpack D, H, z, σ = setup

    rB = sqrt(σ^2 + (z - D - H)^2 )
    rT = sqrt(σ^2 + (z - D)^2     )
    
    r1 = (D < z && z < D+H) ? σ : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)
    
    LineKernelParams(r1=r1, r2=r2, r3=r3, ω=ω, σ=σ, ϵ=ϵ)
end

I_stp(s, D, H, z) = erf(s * (z-D)) - erf(s * (z-D-H))
"""
Compute the number of steps N that can be skipped for a line to point
"""
function compute_N_line(ϵ, Nt, setup::SegmentToPoint, params::Constants) 
    @unpack σ, D, H, z = setup
    @unpack α, kg, Δt = params
    
    f = let kg=kg, σ=σ, α=α, Δt=Δt, D=D, H=H, z=z, ϵ=ϵ
        N -> 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * I_stp(s, D, H, z), 1/sqrt(4*α*Δt*N), Inf)[1] - ϵ
    end

    if f(1)*f(Nt) > 0 return Nt end
    problem = ZeroProblem(f, 10σ^2)
    sol = solve(problem, xtol=1.)
    Int(floor(sol))
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points_line!(ζ, W, N, No, ϵ, ϵ´, n, presetup::SegmentToPoint, params::Constants)
    @unpack Δt, α, rb, kg, Δt̃  = params
    @unpack D, H, z, σ = presetup

    z_int = log((z-D + sqrt(σ^2 + (z-D)^2)) / (z-D-H + sqrt(σ^2 + (z-D-H)^2))) 
    a = 0.
    f = let z_int=z_int, ϵ=ϵ, ϵ´=ϵ´, N=N, No=No, Δt̃=Δt̃
        b -> (gamma(0, b^2*(N-1)*Δt̃) - gamma(0, b^2*(No-1)*Δt̃))/2 * (z_int - ϵ´) + ϵ´ * log((No-1)/(N-1)) - 4π^2*kg * ϵ
    end
    problem = ZeroProblem(f, sqrt(-log(ϵ) / (N*Δt̃)))

    b = abs(solve(problem))
    if b == 0 || isnan(b) b = sqrt(-log(ϵ) / (N*Δt̃)) end
    params = LineKernelParams(presetup, 0., ϵ)

    guide = let N=N, No=No, Δt̃=Δt̃, rb=rb, params=params
        @unpack r1, r2, r3 = params
        ω1 = r1
        ω2 = r1 == r2 ? r3 : r2
        aa(ζ) = sin(ω1/rb * ζ) / ζ / ω1 
        bb(ζ) = sin(ω2/rb * ζ) / ζ / ω2
        ζ -> (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * (aa(ζ)-bb(ζ)) * (1 - exp(-ζ^2*Δt̃)) / ζ
    end
    _, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    sort!(segbuf, by=x->x.a)
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
    return Nζ
end

function compute_kernel_line(params::LineKernelParams; containers)
    @unpack r1, r2, r3, σ, ω, ϵ = params

    I = 0.
    H = let ω=ω, σ=σ
        t -> sin(ω * σ * cosh(t))
    end

    if r2 != r3 
        a = find_integration_interval(ϵ/3, r2, r3, σ, ω, containers)
        I1 = asymptotic(a, r3, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r2/σ), acosh(a/σ), atol=ϵ/3)
        I += I1 + I2
    end

    if r1 != r2
        a = find_integration_interval(ϵ/6, r1, r2, σ, ω, containers)
        I1 = asymptotic(a, r2, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r1/σ), acosh(a/σ), atol=ϵ/6)
        I += 2 * (I1 + I2)
    end
    return I 
end

#=
function bakhalov_discretization(N, setup, params::Constants)
    @unpack D, H, z = setup
    @unpack Δt̃, rb, kg = params
    σ = rb
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
=#

function minimum_distance(source, target)
    σ = compute_distance_2D(source, target)
    if target.D > source.D + source.H 
        dist = sqrt(σ^2 + (target.D - source.D - source.H)^2)
    elseif target.D + target.H < source.D 
        dist = sqrt(σ^2 + (target.D + target.H - source.D )^2)
    else 
        dist = σ
    end
    return dist
end

function minimum_distance_line_to_point(source, target)
    σ = compute_distance_2D(source, target)
    z = target.D + target.H / 2
    if source.D > z
        dist = sqrt(σ^2 + (source.D - z)^2)
    elseif source.D + source.H < z
        dist = sqrt(σ^2 + (source.D + source.H - z)^2)
    else 
        dist = σ
    end
    return dist
end

function get_representative_ltp(sources)
    min_dist = Inf
    setup = SegmentToPoint(D=0., H=0., z=0., σ=0.)
    for (i, s) in enumerate(sources)
        for t in @views sources[1:i-1]
            dist = min(minimum_distance_line_to_point(s, t), minimum_distance_line_to_point(t, s))
            if dist < min_dist
                min_dist = dist
                setup = SegmentToPoint(D=s.D, H=s.H, z=t.D+t.H/2, σ=compute_distance_2D(s, t))
            end
        end
    end
    return setup, min_dist
end

function choose_blocks(::SegmentToPoint, sources, Nt, ϵ, constants)
    @unpack kg, α, Δt, Δt̃, rb = constants
    ratio = 3.
    setup, min_dist = get_representative_ltp(sources)
    @unpack σ, D, H, z = setup

    dline = let σ=σ, D=D, H=H, z=z
        (h, N) -> begin
            if N <= 0 return 0. end
            if D <= z <= D+H
                return 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * abs(2 * erf(s*h) + erf(s*(z-D-H)) - erf(s*(z-D))), 1/sqrt(4*α*Δt*N), Inf, atol=ϵ/100)[1]
            elseif z < D
                return 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * abs(erf(s*(z-D-h)) - erf(s*(z-D-H))), 1/sqrt(4*α*Δt*N), Inf, atol=ϵ/100)[1]
            else 
                return 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * abs(erf(s*(z-D-H+h)) - erf(s*(z-D))), 1/sqrt(4*α*Δt*N), Inf, atol=ϵ/100)[1]
            end
        end
    end

    Nmin = compute_N(min_dist, ϵ, constants)
    N = [Int(floor(find_zero(N -> dline(0, N) - ϵ, isinf(Nmin) ? Nt : Nmin)))]
    ND = Float64[min_dist]
    i = 1

    offset = 0.
    if z < D  
        offset =  D - z
    elseif z > D + H
        offset = z - D - H
    end

    while true
        d = sqrt(σ^2 + offset^2) * ratio^i
        h = sqrt(d^2 - σ^2) - offset
        if (D <= z <= D+H && h < H/2) || h < H 
            Nd = compute_N(d, ϵ, constants)
            Nr = Int(floor(find_zero(N -> dline(h, N) - ϵ, isinf(Nd) ? Nt : Nd)))
        else 
            Nr = compute_N(d, ϵ, constants)
        end
        i += 1
        if Nr < N[end] continue end
        if Nr > Nt || Nt * 0.75 < Nr < Nt break end
        push!(N, Nr)
        push!(ND, d)
    end
    if N[end] < Nt
        push!(N, Nt)
        edge = max(abs(z-D), abs(z-D-H))
        maxh = find_zero(h -> dline(h, Nt) - ϵ, ND[end])
        push!(ND, sqrt(σ^2 + min(edge^2, maxh^2)))
    end
    return N, ND
end 

function compute_distance(::SegmentToPoint, source, target, params, ϵ, Nt)
    dist = minimum_distance(source, target)
    N_r = compute_N(dist, ϵ, params) 
    dist, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToPoint; sources, N, ND, n, ϵ, ϵ´, constants)
    K = length(N) - 1
    σ = ND[1]
    maxH = maximum(map(e -> e.H, sources))
    for i in 1:K
        H = min(2 * sqrt(ND[i+1]^2 - σ^2), maxH)
        setup = SegmentToPoint(D=0., H=H, z=H/2, σ=σ)
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ, ϵ´, n, setup, constants)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToPoint; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers, ND, ranges)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        block = findfirst(range -> k in range, ranges)
        z_eval = target.D + target.H/2
        if z_eval < source.D  
            h = sqrt(ND[block+1]^2 - σ^2) - (source.D - z_eval)
            D_eval = source.D
            H_eval = h
        elseif z_eval > source.D + source.H
            h = sqrt(ND[block+1]^2 - σ^2) - (z_eval - source.D - source.H)
            D_eval = source.D + source.H - h
            H_eval = h
        else
            h = sqrt(ND[block+1]^2 - σ^2)
            D_eval = max(source.D, z_eval - h)
            H_eval = min(source.D + source.H, z_eval + h) - D_eval
        end
        setup = SegmentToPoint(D=D_eval, H=H_eval, z=z_eval, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ)
        @inbounds HM[k, i, j] = C * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_line(lineparams; containers=containers)
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end
