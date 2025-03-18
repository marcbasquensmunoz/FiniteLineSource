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
function compute_ζ_points_line!(ζ, W, N, No, ϵ, n, presetup::SegmentToPoint, params::Constants, containers::AsymptoticContainers)
    @unpack Δt, α, rb, kg, Δt̃  = params
    @unpack D, H, z, σ = presetup
    
    #=
    heatwave = let D=D, H=H, z=z, kg=kg, N=N, Δt=Δt, α=α, ϵ=ϵ
        σ -> 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * I_stp(s, D, H, z), 1/sqrt(4*α*Δt*N), Inf)[1] - ϵ
    end
    problem = ZeroProblem(heatwave, sqrt(N))
    sol = solve(problem, xatol=1e-2)
    if isnan(sol) sol = rb end
    =#

    z_int = log((z-D + sqrt(σ^2 + (z-D)^2)) / (z-D-H + sqrt(σ^2 + (z-D-H)^2)))

    a = 0.
    f = let z_int=z_int, ϵ=ϵ, N=N, Δt̃=Δt̃
        b -> 2ϵ/z_int - (gamma(0, b^2*N*Δt̃) - gamma(0, b^2*(N+1)*Δt̃))
    end
    problem = ZeroProblem(f, sqrt(-log(ϵ) / (N*Δt̃)))
    #@show solve(problem), sqrt(-log(ϵ) / (N*Δt̃))
    b = 1.03 * solve(problem)
    if b == 0 || isnan(b) b = sqrt(-log(ϵ) / (N*Δt̃)) end
    
    #setup = SegmentToPoint(D=presetup.D, H=presetup.H, z=presetup.z, σ=σ)
    params = LineKernelParams(presetup, 0., ϵ)
    #@show b, N, No

    guide = let N=N, No=No, Δt̃=Δt̃, rb=rb, params=params, containers=containers, setup=presetup, ϵ=ϵ
        @unpack r1, r2, r3 = params
        ω1 = r1
        ω2 = r1 == r2 ? r3 : r2
        aa(ζ) = sin(ω1/rb * ζ) / ζ / ω1 
        bb(ζ) = sin(ω2/rb * ζ) / ζ / ω2
        ζ -> (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * (-bb(ζ)+aa(ζ)) * (1 - exp(-ζ^2*Δt̃)) / ζ
        #=
        ζ -> begin  
            int_sin = compute_kernel_line(LineKernelParams(setup, ζ/rb, ϵ), containers=containers)
            return (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * int_sin * (1 - exp(-ζ^2*Δt̃)) / ζ
        end 
        =#
    end
    @show b
    I, E, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    #@show length(segbuf)
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
        a = find_integration_interval(ϵ/3, r1, r2, σ, ω, containers)
        I1 = asymptotic(a, r2, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r1/σ), acosh(a/σ), atol=ϵ/3)
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
function choose_blocks(setup::SegmentToPoint, distances, Nr, Nt, ϵ, constants)
    @unpack kg, α, Δt = constants
    σ = minimum(filter(e -> e != 0, distances))
    H = 150.
    D = 0.
    z = 75.
    ratio = 2.
    max_dist = sqrt(H^2 + σ^2)
    max_N = compute_N(max_dist, ϵ, constants)
    N = Int[compute_N(σ, ϵ, constants)]
    ND = Float64[σ]
    i = 1
    while true
        d = σ * ratio^i
        @show d


        h = sqrt(d^2 - σ^2)
        if h < H/2
            #L = H/2 - d
            #dline(L, N) = -1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * (erf(s*(z-D-L)) - erf(s*(z-D-H+L)) - erf(s*(z-D)) + erf(s*(z-D-H))), 1/sqrt(4*α*Δt*N), Inf)[1]
            dline(h, N) = N <= 0 ? 0 : 1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * abs(erf(s*(z-D-H/2+h)) - erf(s*(z-D-H/2-h)) - erf(s*(z-D)) + erf(s*(z-D-H))), 1/sqrt(4*α*Δt*N), Inf, atol=ϵ)[1]
            Nr = Int(floor(find_zero(N -> dline(h, N) - ϵ, compute_N(d, ϵ, constants))))
            @show Nr
        else 
            Nr = compute_N(d, ϵ, constants)
        end
        #Nr = 868
        #Nr = compute_N(d, ϵ, constants)
        i += 1
        if Nr < N[end] continue end
        if Nr > Nt || Nt * 0.75 < Nr < Nt break end
        #@show Nr, d
        push!(N, Nr)
        push!(ND, d)
    end
    #push!(N, max_N)
    #push!(ND, max_dist)
    if N[end] < Nt
        push!(N, Nt)
        push!(ND, sqrt(H^2/4 + σ^2))
    end
    #@show N, ND
    return N, ND
end 

function compute_distance(::SegmentToPoint, source, target, params, ϵ, Nt)
    σ = compute_distance_2D(source, target)
    #setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H / 2 , σ=σ)
    #N_r = compute_N_line(ϵ, Nt, setup, params)
    if target.D > source.D + source.H 
        dist = sqrt(σ^2 + (target.D - source.D - source.H)^2)
    elseif target.D + target.H < source.D 
        dist = sqrt(σ^2 + (target.D + target.H - source.D )^2)
    else 
        dist = σ
    end

    N_r = compute_N(dist, ϵ, params) 
    dist, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToPoint; sources, N, ND, n, ϵ, constants, containers)
    K = length(N) - 1
    σ = ND[1]
    for i in 1:K
        H = min(2 * sqrt(ND[i+1]^2 - σ^2), maximum(map(e -> e.H, sources)))
        setup = SegmentToPoint(D=0., H=H, z=H/2, σ=σ)
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ/(length(N)-1), n, setup, constants, containers)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToPoint; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers, N, ND, ranges)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        block = findfirst(range -> k in range, ranges)
        h = sqrt(ND[block+1]^2 - σ^2)
        z_eval = target.D + target.H/2
        D_eval = max(source.D, z_eval - h)
        H_eval = min(source.D + source.H, z_eval + h) - D_eval
        setup = SegmentToPoint(D=D_eval, H=H_eval, z=z_eval, σ=σ)
        #setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H/2, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ)
        @inbounds HM[k, i, j] = C * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_line(lineparams; containers=containers)
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end
