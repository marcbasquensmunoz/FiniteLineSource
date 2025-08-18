@with_kw struct LineToLineIntegralParams{T <: Number} @deftype T
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

function LineToLineIntegralParams(setup::SegmentToSegment, ω, ϵ, h=nothing)
    @unpack D1, H1, D2, H2, σ = setup

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : σ
    r2 = H2 > H1 ? (D2 > D1 ? rLL : σ) : (D2 + H2 > D1 + H1 ? rUR : σ)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : σ) : (D2 > D1 ? rLL : σ)
    r4 = D2 + H2 > D1 ? rUL : σ

    α1 = abs(D1 - D2 + H1)
    α2 = min(H1, H2)
    α3 = abs(D2 - D1 + H2)

    d = isnothing(h) ? Inf : sqrt(σ^2 + h^2)
    LineToLineIntegralParams(r1=min(r1, d), r2=min(r2, d), r3=min(r3, d), r4=min(r4, d), α1=α1, α2=α2, α3=α3, ω=ω, σ=σ, ϵ=ϵ)
end

function compute_distance(::SegmentToSegment, source, target, constants::Constants, ϵ, Nt, Q)
    dist = minimum_distance(source, target)
    new_setup = SegmentToSegment(D1=source.D, H1=source.H, D2=target.D, H2=target.H, σ=dist)
    N_r = compute_N_for_line_to_line_range(0., new_setup, constants, ϵ, Nt, Nt, Q)
    dist, N_r
end

function constant_integral(setup::SegmentToSegment, constants::Constants, N) 
    @unpack D1, H1, D2, H2, σ = setup
    @unpack Δt̃, kg = constants
    rb = σ
    r(z1, z2) = sqrt(rb^2 + (z1 - z2)^2)
    quadgk(z1 -> quadgk(z2 -> erf(r(z1, z2)/rb/sqrt(4*N*Δt̃)) / r(z1, z2), D2, D2+H2)[1], D1, D1+H1)[1] / (4π * kg * H2)
end

self_setup(::SegmentToSegment, source) = SegmentToSegmentOld(D1=source.D, H1=source.H, D2=source.D, H2=source.H, σ=source.rb)
image(s::SegmentToSegmentOld) = SegmentToSegmentOld(D1=-s.D1-s.H1, H1=s.H1, D2=s.D2, H2=s.H2, σ=s.σ)

function get_representative_ltl(sources, image_strength)
    min_dist = Inf
    setup = SegmentToSegment(D1=0., H1=0., D2=0., H2=0., σ=0.)
    for (i, s) in enumerate(sources)
        for t in @views sources[1:i-1]
            dist = minimum_distance(s, t)
            if dist < min_dist
                min_dist = dist
                setup = SegmentToSegment(D1=s.D, H1=s.H, D2=t.D, H2=t.H, σ=compute_distance_2D(s, t), image_strength = image_strength)
            end
        end
    end
    return setup, min_dist
end

Θ(x) = x >= 0 ? 1. : 0.
int_exp_ierf(x, L, σ) = x/2 * gamma(0, σ^2*L^2) - exp(-σ^2*L^2) / L / sqrt(π) + σ*erfc(σ*L) + exp(-(σ^2+x^2)*L^2) / L / sqrt(π) - sqrt(σ^2+x^2)*erfc(sqrt(σ^2+x^2)*L) 

function I_L2L(h, N, setup, constants, ϵ)
    if N <= 0 return 0. end
    @unpack α, kg, Δt = constants
    @unpack D1, H1, D2, H2, σ, image_strength = setup
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
    real = quadgk(integrand, 1/sqrt(4α*Δt*N), Inf, atol=ϵ)[1]

    error = real - IA

    # Image
    if image_strength != 0.
        I_h_image(s) = (
            1/s * (max(h-D1-H1, min(D2+H2, h-D1)) - max(h-D1-H1, min(D2, h-D1))) * erf(s*h)
            + 1/s^2 * ( ierf(s*(min(h-D1-H1+D2+H2, 2*(D2+H2)))) - ierf(s*(min(h-D1-H1+D2+H2, 2D2+H2))) )
            + 1/s^2 * ( ierf(s*(min(D2+D1, h))) - ierf(s*(min(D2+H2+D1, h))) )
        )
        I_full_image(s) = (ierf(s*(D2+H2+D1+H1)) + ierf(s*(D1+D2)) - ierf(s*(D2+H2+D1)) - ierf(s*(D2+D1+H1))) / s^2

        image_error = quadgk(s -> exp(-s^2*σ^2) * (I_full_image(s)-I_h_image(s)), 1/sqrt(4α*Δt*N), Inf, atol=ϵ)[1]
        error += image_strength * image_error
    end
    
    return 1 / (4π*kg * H2) * abs(error)
end

function compute_N_for_line_to_line_range(h, setup, constants, ϵ, Nt, Ndef, Q)
    try 
        return Int(floor(find_zero(N -> I_L2L(h, N, setup, constants, ϵ) - ϵ/Q, isinf(Ndef) ? Nt : Ndef)))
    catch 
        return Ndef
    end
end

function choose_blocks(setup::SegmentToSegment, sources, Nt, ϵ, constants, Q)
    if length(sources) == 1
        return [Nt], [sources[1].rb * 100]
    end
    @unpack kg, α, Δt, Δt̃ = constants
    @unpack D1, H1, D2, H2, σ, image_strength = setup

    ratio = 3.
    setup, min_dist = get_representative_ltl(sources, image_strength)

    Nmin = compute_N(min_dist, ϵ, constants, Nt, Q)

    if I_L2L(0., Nt, setup, constants, ϵ) < ϵ/Q
        N = Int[Nt]
    else 
        N = [compute_N_for_line_to_line_range(0., setup, constants, ϵ, Nt, Nmin, Q)]
    end
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
        if h < H1+H2
            if I_L2L(h, Nt, setup, constants, ϵ) < ϵ/Q
                break
            end
            Nd = compute_N(d, ϵ, constants, Q)
            Nr = compute_N_for_line_to_line_range(h, setup, constants, ϵ, Nt, Nd, Q)
        else 
            Nr = compute_N(d, ϵ, constants, Q)
        end
        i += 1
        if Nr < N[end] continue end
        if Nr > Nt || Nt * 0.75 < Nr < Nt break end
        if !(Nr in N) && N[end] < Nr * 0.75
            push!(N, Nr)
            push!(ND, d)
        end
    end
    if N[end] < Nt
        push!(N, Nt)
        edge = max(abs(D1+H1-D2), abs(D2+H2-D1))
        local maxh
        try
            maxh = find_zero(h -> I_L2L(h, Nt, setup, constants, ϵ) - ϵ/Q, ND[end])
        catch
            maxh = edge
        end
        push!(ND, sqrt(σ^2 + min(edge^2, maxh^2)))
    end
    
    zero_indices = findall(==(0), N)
    deleteat!(N, zero_indices)
    deleteat!(ND, zero_indices)
    #=
    push!(N, 1)
    push!(ND, compute_r(1, ϵ, constants, Q))
    perm = sortperm(N)
    return N[perm], ND[perm]
    =#
    return N, ND
end 

zero_freq(a, b, σ) = log( ( (b+sqrt(b^2-σ^2)) * (-a+sqrt(a^2-σ^2)) ) / ( (a+sqrt(a^2-σ^2)) * (-b+sqrt(b^2-σ^2)) ) ) / 2

function compute_double_line_integral_part(params::LineToLineIntegralParams; containers)
    @unpack r1, r2, r3, r4, ω, ϵ, α1, α2, α3, σ = params

    if ω == 0.
        return r2 + r3 - r1 - r4 + zero_freq(r1, r2, σ) * α1 + zero_freq(r2, r3, σ) * α2 + zero_freq(r3, r4, σ) * α3
    end
        
    I = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω

    H = let ω=ω, σ=σ
        t -> sin(ω * σ * cosh(t))
    end

    #@show r1, r2, r3, r4
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

function compute_double_line_integral(params::LineToLineIntegralParams, paramsT::LineToLineIntegralParams; containers)
    I1 = compute_double_line_integral_part(params, containers=containers)
    I2 = compute_double_line_integral_part(paramsT, containers=containers)
    (I1+I2)
end

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points_line_to_line!(ζ, W, N, No, ϵ, ϵ´, n, presetup::SegmentToSegment, constants::Constants, Q)
    @unpack Δt, α, rb, kg, Δt̃ = constants
    @unpack D1, H1, D2, H2, σ = presetup
    β = let σ=σ
        d -> sqrt(σ^2+d^2) + d * log(sqrt(σ^2+d^2) - d)
    end
    z_int = β(D1+H1-D2-H2) + β(D1-D2) - β(D1-H2-D2) - β(D1+H1-D2)

    a = 0.
    if N != 0
        f = let z_int=z_int, ϵ=ϵ, ϵ´=ϵ´, N=N, No=No, Δt̃=Δt̃, kg=kg
            b -> (gamma(0, b^2*(N-1)*Δt̃) - gamma(0, b^2*(No-1)*Δt̃)) * (z_int - ϵ´) + ϵ´ * log((No-1)/(N-1)) - 4*π^2*kg * ϵ/Q
        end    
        b0 = sqrt(-log(ϵ/Q) / (N*Δt̃))
        problem = ZeroProblem(f, b0)
        sol_b = abs(solve(problem))
        b = isnan(sol_b) || sol_b == 0 || sol_b > 10 ? b0 : sol_b
    else 
        b = 10.
    end

    params = LineToLineIntegralParams(presetup, 0., ϵ/Q)
    paramsT = LineToLineIntegralParams(transpose(presetup), 0., ϵ/Q)
    ω1 = min(params.r1, paramsT.r1)
    ω2 = max(params.r4, paramsT.r4)

    guide = let N=N, No=No, ω1=ω1, ω2=ω2, Δt̃=Δt̃, kg=kg, rb=rb, Q=Q
        ζ -> Q/(2 * π^2 * kg) * (acosh(ω2/σ) - acosh(ω1/σ)) * (No-N) * (exp(-ζ^2*N*Δt̃) - exp(-ζ^2*No*Δt̃)) * (sin(ζ*ω2/rb) - sin(ζ*ω1/rb)) #=* (1 - exp(-ζ^2*Δt̃))=# / ζ
    end
    _, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=Q*ϵ)
    
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

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToSegment; sources, N, ND, n, ϵ, ϵ´, constants, Q=1.)
    K = length(N) - 1
    σ = ND[1]
    for i in 1:K
        H = min(2 * sqrt(ND[i+1]^2 - σ^2), maximum(map(e -> e.H, sources)))
        setup = SegmentToSegment(D1=0., H1=H, D2=0., H2=H, σ=σ)
        N_block = compute_ζ_points_line_to_line!(ζ, W, N[i], N[i+1], ϵ, ϵ´, n, setup, constants, Q)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, setup::SegmentToSegment; ζ, W, expt, sources, distances, constants::Constants, ϵ, containers, ND, ranges)
    @unpack kg = constants
    @unpack image_strength = setup

    C = 1 / (2π^2*kg)
    for j in eachindex(sources), i in 1:j, (k, ζζ) in enumerate(ζ)
        rb = sources[j].rb
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        block = findfirst(range -> k in range, ranges)
        
        if σ > ND[block+1] continue end
        h = sqrt(ND[block+1]^2 - σ^2)

        setup = SegmentToSegment(D1=source.D, H1=source.H, D2=target.D, H2=target.H, σ=σ)

        lineparams = FiniteLineSource.LineToLineIntegralParams(setup, ζζ/rb, ϵ, h)
        lineparamsT = FiniteLineSource.LineToLineIntegralParams(transpose(setup), ζζ/rb, ϵ, h)

        interaction = compute_double_line_integral(lineparams, lineparamsT; containers=containers) 

        if image_strength != 0.
            setup_image = SegmentToSegment(D1=-source.D-source.H, H1=source.H, D2=target.D, H2=target.H, σ=σ)

            lineparams_image = FiniteLineSource.LineToLineIntegralParams(setup_image, ζζ/rb, ϵ, h)
            lineparamsT_image = FiniteLineSource.LineToLineIntegralParams(transpose(setup_image), ζζ/rb, ϵ, h)

            image_interaction = compute_double_line_integral(lineparams_image, lineparamsT_image; containers=containers)
            interaction += image_strength * image_interaction
        end

        @inbounds HM[k, i, j] = C / target.H * W[k] * (1 - expt[k]) / ζ[k] * interaction
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end