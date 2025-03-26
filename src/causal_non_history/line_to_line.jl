
@with_kw struct LineToLineKernelParams{T <: Number} @deftype T
    r1
    r2
    r3
    r4
    α1
    α2
    α3
    ω
    σ
    ϵ
end

function LineToLineKernelParams(setup::SegmentToSegment, ω, ϵ, h=nothing)
    @unpack D1, H1, D2, H2, σ = setup

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : σ
    r2 = H2 > H1 ? (D2 > D1 ? rLL : σ) : (D2 + H2 > D1 + H1 ? rUR : σ)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : σ) : (D2 > D1 ? rLL : σ)
    r4 = D2 + H2 > D1 ? rUL : σ

    α1 = D1 - D2 + H1
    α2 = min(H1, H2)
    α3 = D2 - D1 + H2

    d = isnothing(h) ? Inf : sqrt(σ^2 + h^2)
    LineToLineKernelParams(r1=min(r1, d), r2=min(r2, d), r3=min(r3, d), r4=min(r4, d), α1=α1, α2=α2, α3=α3, ω=ω, σ=σ, ϵ=ϵ)
end

function minimum_distance_line_to_line(source, target)
    σ = compute_distance_2D(source, target)
    if target.D > source.D + source.H
        dist = sqrt(σ^2 + (target.D - source.D - source.H)^2)
    elseif source.D > target.D + target.H
        dist = sqrt(σ^2 + (source.D - target.D - target.H)^2)
    else 
        dist = σ
    end
    return dist
end

function get_representative_ltl(sources)
    min_dist = Inf
    setup = SegmentToSegment(D1=0., H1=0., D2=0., H2=0., σ=0.)
    for (i, s) in enumerate(sources)
        for t in @views sources[1:i-1]
            dist = minimum_distance_line_to_line(s, t)
            if dist < min_dist
                min_dist = dist
                setup = SegmentToSegment(D1=s.D, H1=s.H, D2=t.D, H2=t.H, σ=compute_distance_2D(s, t))
            end
        end
    end
    return setup, min_dist
end


int_exp_ierf(x, L, σ) = x/2 * gamma(0, σ^2*L^2) - exp(-σ^2*L^2) / L / sqrt(π) + σ*erfc(σ*L) + exp(-(σ^2+x^2)*L^2) / L / sqrt(π) - sqrt(σ^2+x^2)*erfc(sqrt(σ^2+x^2)*L) 

ierf(x) = x*erf(x) - 1/sqrt(π)*(1- exp(-x^2))
function I_L2L(h, N, setup, constants)
    if N <= 0 return 0. end
    @unpack α, kg, Δt = constants
    @unpack D1, H1, D2, H2, σ = setup
    P = min(H2, D2 + H2 - D1 - h) + min(H2, D1 + H1 - D2 - h)
    M1 = min(D2 + H2 - D1, h)
    M2 = min(D1 + H1 - D2, h)
    A1 = D2 + H2 - D1
    A2 = D1 + H1 - D2

    threshold = 10*sqrt(4α*Δt*N)
    fA = let threshold=threshold, A1=A1, A2=A2
        s -> begin
            res = 0.
            if A1 < threshold 
                res += ierf(s * A1)
            end
            if A2 < threshold
                res += ierf(s * A2)
            end 
            res  
        end
    end
    IA = 0.
    if A1 > threshold 
        IA += int_exp_ierf(A1, 1/sqrt(4α*Δt*N), σ)
    end
    if A2 > threshold
        IA += int_exp_ierf(A2, 1/sqrt(4α*Δt*N), σ) 
    end    

    integrand(s) = exp(-σ^2 * s^2) / s * ((ierf(s * M1) + ierf(s * M2) - fA(s)) / s + erf(s * h) * P)      
    return 1 / (4π*kg * H2) * abs(quadgk(integrand, 1/sqrt(4α*Δt*N), Inf)[1] - IA)
end

function choose_blocks(::SegmentToSegment, sources, Nt, ϵ, constants)
    @unpack kg, α, Δt, Δt̃, rb = constants

    ratio = 3.
    setup, min_dist = get_representative_ltl(sources)
    @unpack D1, H1, D2, H2, σ = setup

    N = [Int(floor(find_zero(N -> I_L2L(0., N, setup, constants) - ϵ, compute_N(min_dist, ϵ, constants))))]
    ND = Float64[min_dist]
    i = 1

    offset = 0.
    if D1 + H1 < D2  
        offset =  D2 - D1 - H1
    elseif D1 > D2 + H2
        offset = D1 - D2 - H2
    end

    while true
        d = sqrt(σ^2 + offset^2) * ratio^i
        h = sqrt(d^2 - σ^2) - offset
        if h < H2/2
            Nr = Int(floor(find_zero(N -> I_L2L(h, N, setup, constants) - ϵ, compute_N(d, ϵ, constants), xatol=1.)))
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
        push!(ND, sqrt(σ^2 + max((D1-D2-H2)^2, (D2-D1-H1)^2)))
    end
    return N, ND
end 

zero_freq(a, b, σ) = log( ( (b+sqrt(b^2-σ^2)) * (-a+sqrt(a^2-σ^2)) ) / ( (a+sqrt(a^2-σ^2)) * (-b+sqrt(b^2-σ^2)) ) ) / 2

function compute_kernel_double_line_part(params::LineToLineKernelParams; containers)
    @unpack r1, r2, r3, r4, ω, ϵ, α1, α2, α3, σ = params

    if ω == 0.
        return r2 + r3 - r1 - r4 + zero_freq(r1, r2, σ) * α1 + zero_freq(r2, r3, σ) * α2 + zero_freq(r3, r4, σ) * α3
    end
        
    I = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω

    H = let ω=ω, σ=σ
        t -> sin(ω * σ * cosh(t))
    end

    if r3 != r4 
        I3 = 0.
        if r1 == r3
            a = find_integration_interval(ϵ/(3*α3), r3, r4, σ, ω, containers)
            I3 += asymptotic(a, r4, σ, ω, containers)
        else 
            a = r4
        end
        I3 += quadgk(H, acosh(r3/σ), acosh(a/σ), atol=ϵ/(3*α3))[1]
        I += α3 * I3
    end
    if r2 != r3 
        I2 = 0.
        if r1 == r2
            a = find_integration_interval(ϵ/(3*α2), r2, r3, σ, ω, containers)
            I2 += asymptotic(a, r3, σ, ω, containers)
        else
            a = r3
        end
        I2 += quadgk(H, acosh(r2/σ), acosh(a/σ), atol=ϵ/(3*α2))[1]
        I += α2 * I2
    end
    if r1 != r2
        a = find_integration_interval(ϵ/(3*α1), r1, r2, σ, ω, containers)
        I1 = asymptotic(a, r2, σ, ω, containers)
        I2, _ = quadgk(H, acosh(r1/σ), acosh(a/σ), atol=ϵ/(3*α1))
        I += α1 * (I1 + I2)
    end
    return I
end

function compute_kernel_double_line(params::LineToLineKernelParams, paramsT::LineToLineKernelParams; containers)
    I1 = compute_kernel_double_line_part(params, containers=containers)
    I2 = compute_kernel_double_line_part(paramsT, containers=containers)
    (I1+I2)
end

function compute_distance(::SegmentToSegment, source, target, params, ϵ, Nt)
    dist = minimum_distance(source, target)
    N_r = compute_N(dist, ϵ, params) 
    dist, N_r
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points_line_to_line!(ζ, W, N, No, ϵ, ϵ´, n, presetup::SegmentToSegment, constants::Constants)
    @unpack Δt, α, rb, kg, Δt̃ = constants
    @unpack D1, H1, D2, H2, σ = presetup
    β = let σ=σ
        d -> sqrt(σ^2+d^2) + d * log(sqrt(σ^2+d^2) - d)
    end
    z_int = 1/(4π*kg*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1-H2-D2) - β(D1+H1-D2))

    a = 0.
    f = let z_int=z_int, ϵ=ϵ, ϵ´=ϵ´, N=N, No=No, Δt̃=Δt̃
        b -> z_int * (gamma(0, b^2*(N-1)*Δt̃) - gamma(0, b^2*(No-1)*Δt̃))/2 + ϵ´ * sqrt(π)/2 * ( erf(b*sqrt((N-1)*Δt̃)) / sqrt((N-1)*Δt̃) - erf(b*sqrt((No-1)*Δt̃)) / sqrt((No-1)*Δt̃) ) - ϵ
    end    
    problem = ZeroProblem(f, sqrt(-log(ϵ) / (N*Δt̃)))
    sol_b = solve(problem)
    b = isnan(sol_b) || sol_b == 0 ? sqrt(-log(ϵ) / (N*Δt̃)) : sol_b
  
    params = LineToLineKernelParams(presetup, 0., ϵ)
    paramsT = LineToLineKernelParams(transpose(presetup), 0., ϵ)
    ω1 = min(params.r1, paramsT.r1)
    ω2 = max(params.r4, paramsT.r4)

    guide = let N=N, No=No, Δt̃=Δt̃, rb=rb, ω1=ω1, ω2=ω2
        aa(ζ) = sin(ω1/rb * ζ) / ζ / ω1 
        bb(ζ) = sin(ω2/rb * ζ) / ζ / ω2
        ζ -> ζ == 0 ? 0. : (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * (aa(ζ)-bb(ζ)) * (1 - exp(-ζ^2*Δt̃)) / ζ
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

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToSegment; sources, N, ND, n, ϵ, ϵ´, constants)
    K = length(N) - 1
    σ = ND[1]
    for i in 1:K
        H = min(2 * sqrt(ND[i+1]^2 - σ^2), maximum(map(e -> e.H, sources)))
        setup = SegmentToSegment(D1=0., H1=H, D2=0., H2=H, σ=σ)
        N_block = compute_ζ_points_line_to_line!(ζ, W, N[i], N[i+1], ϵ, ϵ´, n, setup, constants)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToSegment; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers, ND, ranges)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)
    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        block = findfirst(range -> k in range, ranges)
        h = sqrt(ND[block+1]^2 - σ^2)
        setup = SegmentToSegment(D1=source.D, H1=source.H, D2=target.D, H2=target.H, σ=σ)
        lineparams = FiniteLineSource.LineToLineKernelParams(setup, ζζ/rb, ϵ, h)
        lineparamsT = FiniteLineSource.LineToLineKernelParams(transpose(setup), ζζ/rb, ϵ, h)

        @inbounds HM[k, i, j] = C / target.H * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_double_line(lineparams, lineparamsT; containers=containers)
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end