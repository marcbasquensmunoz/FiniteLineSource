
@with_kw struct LineToLineKernelParams{T <: Number} @deftype T
    r1
    r2
    r3
    r4
    rs
    α1
    α2
    α3
    ω
    σ
    ϵ
    n1::Int
    n2::Int
    n3::Int
end

function LineToLineKernelParams(setup::SegmentToSegment, ω, ϵ, nmodel)
    @unpack D1, H1, D2, H2, σ = setup

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : σ
    r2 = H2 > H1 ? (D2 > D1 ? rLL : σ) : (D2 + H2 > D1 + H1 ? rUR : σ)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : σ) : (D2 > D1 ? rLL : σ)
    r4 = D2 + H2 > D1 ? rUL : σ

    #rs = r1 + 4π/ω
    rt = r1 == r2 ? (r1 == r3 ? r4 : r3) : r2
    rs = min(r1 + max(2, (rt-r1) * 0.01), rt)

    α1 = D1 - D2 + H1
    α2 = min(H1, H2)
    α3 = D2 - D1 + H2
    n1, n2, n3 = 0, 0, 0

    if r1 != r2
        n1 = N_bound(ϵ / (3α1 * sqrt(r2-rs)), rs, r2, σ, nmodel)
    end

    if r2 != r3
        rl = max(r2, rs)
        n2 = N_bound(ϵ / (3α2 * sqrt(r3-rl)), rl, r3, σ, nmodel)
    end

    if r3 != r4
        rl = max(r3, rs)
        n3 = N_bound(ϵ / (3α3 * sqrt(r4-rl)), rl, r4, σ, nmodel)
    end

    LineToLineKernelParams(r1=r1, r2=r2, r3=r3, r4=r4, rs=rs, α1=α1, α2=α2, α3=α3, ω=ω, σ=σ, ϵ=ϵ, n1=n1, n2=n2, n3=n3)
end


function compute_kernel_double_line_part(params::LineToLineKernelParams; BV1, BV2, BV3)
    @unpack r1, r2, r3, r4, rs, ω, ϵ, α1, α2, α3, σ = params

    h(r) = r < r2 ? α1 : (r < r3 ? α2 : α3)
    I_div = 0.
    I_osc = 0.

    if r1 == σ
        mult = r1 == r2 ? (r2 == r3 ? α3 : α2) : α1
        C = mult * sin(σ*ω)
        h_reg(r) = r == σ ? 0. : (sin(r*ω) * h(r) - C) / sqrt(r^2-σ^2)
        I1, _ = quadgk(h_reg, r1, rs, atol=ϵ/3)
        I2 = C * acoth(rs/sqrt(rs^2-r1^2))
        I_div = I1 + I2
    end

    if r1 != r2
        x1 = BV1.x
        X1 = BV1.X
        P1 = BV1.P
        f1 = BV1.f
        n1 = BV1.n

        m1 = (r2-rs)/2
        c1 = (r2+rs)/2
        @. X1 = m1*x1 + c1

        @. f1 = α1 / sqrt(X1^2-σ^2)
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        I_osc += sqrt(m1*π/(2ω)) * dot(X1, P1, f1)
    end

    if r2 != r3    
        x2 = BV2.x
        X2 = BV2.X
        P2 = BV2.P
        f2 = BV2.f
        n2 = BV2.n

        rl = max(rs, r2)
        m2 = (r3-rl)/2
        c2 = (r3+rl)/2
        @. X2 = m2*x2 + c2

        @. f2 = α2 / sqrt(X2^2-σ^2)
        besselj!(X2, 1/2:(n2+1/2), m2*ω)
        @. X2 = X2 * imag(exp(im*ω*c2) * im^(0:n2)) * (2(0:n2)+1)
        I_osc += sqrt(m2*π/(2ω)) * dot(X2, P2, f2)
    end

    if r3 != r4
        x3 = BV3.x
        X3 = BV3.X
        P3 = BV3.P
        f3 = BV3.f
        n3 = BV3.n

        rl = max(rs, r3)
        m3 = (r4-rl)/2
        c3 = (r4+rl)/2
        @. X3 = m3*x3 + c3

        @. f3 = α3 / sqrt(X3^2-σ^2)
        besselj!(X3, 1/2:(n3+1/2), m3*ω)
        @. X3 = X3 * imag(exp(im*ω*c3) * im^(0:n3)) * (2(0:n3)+1)
        I_osc += sqrt(m3*π/(2ω)) * dot(X3, P3, f3)
    end
    I_osc_linear = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω
    I_div + I_osc + I_osc_linear
end

function compute_kernel_double_line(params::LineToLineKernelParams, paramsT::LineToLineKernelParams; BV1, BV2, BV3)
    I1 = compute_kernel_double_line_part(params, BV1=BV1, BV2=BV2, BV3=BV3)
    I2 = compute_kernel_double_line_part(paramsT, BV1=BV1, BV2=BV2, BV3=BV3)
    (I1+I2)
end

function compute_distance(::SegmentToSegment, source, target, params, ϵ)
    σ = compute_distance_2D(source, target)
    setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H / 2 , σ=σ)
    N_r = compute_N_line(ϵ, setup, params)
    σ, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToSegment; sources, N, n, ϵ, constants, nmodel)
    K = length(N) - 1
    D_eff = sum([source.D for source in sources])/length(sources)
    H_eff = sum([source.H for source in sources])/length(sources)
    z_eff = D_eff + H_eff/2
    for i in 1:K
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ, n, D_eff, H_eff, z_eff, constants, nmodel)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToSegment; ζ, W, expt, sources, distances, constants::Constants, ϵ, nmodel)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    containers = create_bin_containers(bins)

    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]

        setup = SegmentToSegment(D1=source.D, H1=source.H, D2=target.D, H2=target.H, σ=σ)
        lineparams = FiniteLineSource.LineToLineKernelParams(setup, ζζ/rb, ϵ, nmodel)
        lineparamsT =FiniteLineSource.LineToLineKernelParams(transpose(setup), ζζ/rb, ϵ, nmodel)

        bin1 = FiniteLineSource.get_bin(max(lineparams.n1, lineparamsT.n1), bins)
        bin2 = FiniteLineSource.get_bin(max(lineparams.n2, lineparamsT.n2), bins)
        bin3 = FiniteLineSource.get_bin(max(lineparams.n3, lineparamsT.n3), bins)

        @inbounds HM[k, i, j] = C / target.H * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_double_line(lineparams, lineparamsT; BV1=containers[bin1], BV2=containers[bin2], BV3=containers[bin3])
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end