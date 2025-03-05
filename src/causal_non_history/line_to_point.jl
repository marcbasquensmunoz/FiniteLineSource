
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

"""
Compute the number of steps N that can be skipped for a line to point
"""
function compute_N_line(ϵ, setup::SegmentToPoint, params::Constants) 
    @unpack σ, D, H, z = setup
    @unpack α, kg, Δt = params
    f(N) = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z)^2) / sqrt(4α * Δt * N))/(4*π*kg*sqrt(σ^2 + (zp - z)^2)), D, D+H)[1] - ϵ
    problem = ZeroProblem(f, 10σ^2)
    sol = solve(problem)
    Int(floor(sol))
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points_line!(ζ, W, N, No, ϵ, n, D, H, z, params::Constants, containers)
    @unpack Δt, α, rb, kg, Δt̃  = params

    heatwave(σ) = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z)^2) / sqrt(4α * Δt * N))/(4*π*kg*sqrt(σ^2 + (zp - z)^2)), D, D+H)[1] - ϵ
    problem = ZeroProblem(heatwave, sqrt(N))
    σ = solve(problem)

    z_int = log((z-D + sqrt(σ^2 + (z-D)^2)) / (z-D-H + sqrt(σ^2 + (z-D-H)^2)))

    a = 0.
    f(b) = 2ϵ/z_int - (gamma(0, b^2*N*Δt̃) - gamma(0, b^2*(N+1)*Δt̃))
    problem = ZeroProblem(f, sqrt(-log(ϵ) / (N*Δt̃)))
    b = solve(problem)
    if b == 0
        b = sqrt(-log(ϵ) / (N*Δt̃))
    end

    setup = SegmentToPoint(D=D, H=H, z=z, σ=σ)

    int_sin(ζ) = compute_kernel_line(LineKernelParams(setup, ζ/rb, ϵ), containers=containers)
    guide(ζ) = (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * int_sin(ζ) * (1 - exp(-ζ^2*Δt̃)) / ζ
   #_, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)[3]

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
    #H(t) = sin(ω * σ * cosh(t))

    if r2 != r3 
        a = find_integration_interval(ϵ/4, r2, r3, σ, ω, containers)
        I1 = asymptotic(a, r3, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r2/σ), acosh(a/σ), atol=ϵ/4)
        I += I1 + I2
    end

    if r1 != r2
        a = find_integration_interval(ϵ/4, r1, r2, σ, ω, containers)
        I1 = asymptotic(a, r2, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r1/σ), acosh(a/σ), atol=ϵ/4)
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

function compute_distance(::SegmentToPoint, source, target, params, ϵ)
    σ = compute_distance_2D(source, target)
    setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H / 2 , σ=σ)
    N_r = compute_N_line(ϵ, setup, params)
    σ, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToPoint; sources, N, n, ϵ, constants, containers)
    K = length(N) - 1
    D_eff = sum([source.D for source in sources])/length(sources)
    H_eff = sum([source.H for source in sources])/length(sources)
    z_eff = D_eff + H_eff/2
    for i in 1:K
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ, n, D_eff, H_eff, z_eff, constants, containers)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToPoint; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H/2, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ)
        @inbounds HM[k, i, j] = C * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_line(lineparams; containers=containers)
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end
