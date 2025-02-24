
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


function compute_kernel_double_line_part(params::LineToLineKernelParams; x1=nothing, x2=nothing, x3=nothing, X1=nothing, X2=nothing, X3=nothing, f1=nothing, f2=nothing, f3=nothing, P1=nothing, P2=nothing, P3=nothing)
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
        m1 = (r2-rs)/2
        c1 = (r2+rs)/2
        @. X1 = m1*x1 + c1
        n1 = length(X1)-1

        @. f1 = α1 / sqrt(X1^2-σ^2)
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        I_osc += sqrt(m1*π/(2ω)) * dot(X1, P1, f1)
    end

    if r2 != r3    
        rl = max(rs, r2)
        m2 = (r3-rl)/2
        c2 = (r3+rl)/2
        @. X2 = m2*x2 + c2
        n2 = length(X2)-1

        @. f2 = α2 / sqrt(X2^2-σ^2)
        besselj!(X2, 1/2:(n2+1/2), m2*ω)
        @. X2 = X2 * imag(exp(im*ω*c2) * im^(0:n2)) * (2(0:n2)+1)
        I_osc += sqrt(m2*π/(2ω)) * dot(X2, P2, f2)
    end

    if r3 != r4
        rl = max(rs, r3)
        m3 = (r4-rl)/2
        c3 = (r4+rl)/2
        @. X3 = m3*x3 + c3
        n3 = length(X3)-1

        @. f3 = α3 / sqrt(X3^2-σ^2)
        besselj!(X3, 1/2:(n3+1/2), m3*ω)
        @. X3 = X3 * imag(exp(im*ω*c3) * im^(0:n3)) * (2(0:n3)+1)
        I_osc += sqrt(m3*π/(2ω)) * dot(X3, P3, f3)
    end
    I_osc_linear = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω
    I_div + I_osc + I_osc_linear
end

function compute_kernel_double_line_part(ζ, rb, setup::SegmentToSegment; x, X, P, f, atol)
    @unpack D1, H1, D2, H2, σ = setup
    rmin = σ

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin

    α1 = D1 - D2 + H1
    α2 = min(H1, H2)
    α3 = D2 - D1 + H2

    ω = ζ/rb
    h1(r) = r < r2 ? α1 : (r < r3 ? α2 : α3)
    h2(r) = r < r2 ? 1. : (r < r3 ? 0. : -1.)

    split = r1
    I_div = 0.

    if r1 == σ
        split = r1 + min(π/ω, (r4-r1)*0.1)
        mult = r1 == r2 ? (r2 == r3 ? α3 : α2) : α1
        C = mult * sin(σ*ω)
        h_reg(r) = r == σ ? 0. : (sin(r*ω) * h1(r) - C) / sqrt(r^2-σ^2)
        I1, _ = quadgk(h_reg, r1, split, atol=atol)
        I2 = C * acoth(split/sqrt(split^2-r1^2))
        I_div = I1 + I2
    end

    m = (r4-split)/2
    c = (r4+split)/2

    n = length(x)-1
    @. X = m*x + c
    @. f = h1(X) / sqrt(X^2-σ^2)

    besselj!(X, 1/2:(n+1/2), m*ω)
    @. X = X * imag(exp(im*ω*c) * im^(0:n)) * (2(0:n)+1)

    I_osc = sqrt(m*π/(2ω)) * X' * P * f
    I_osc_linear = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω

    I_div + I_osc + I_osc_linear
end

function compute_kernel_double_line(ζ, rb, setup::SegmentToSegment; x, X, P, f, atol=1e-8)
    I1 = compute_kernel_double_line_part(ζ, rb, setup, x=x, X=X, P=P, f=f, atol=atol)
    tranposed_setup = SegmentToSegment(D1=setup.D2, H1=setup.H2, D2=setup.D1, H2=setup.H1, σ=setup.σ)
    I2 = compute_kernel_double_line_part(ζ, rb, tranposed_setup, x=x, X=X, P=P, f=f, atol=atol)
    (I1+I2)/setup.H2
end

function compute_kernel_double_line(params::LineToLineKernelParams, paramsT::LineToLineKernelParams; x1=nothing, x2=nothing, x3=nothing, X1=nothing, X2=nothing, X3=nothing, f1=nothing, f2=nothing, f3=nothing, P1=nothing, P2=nothing, P3=nothing)
    I1 = compute_kernel_double_line_part(params, x1=x1, x2=x2, x3=x3, X1=X1, X2=X2, X3=X3, f1=f1, f2=f2, f3=f3, P1=P1, P2=P2, P3=P3)
    I2 = compute_kernel_double_line_part(paramsT, x1=x1, x2=x2, x3=x3, X1=X1, X2=X2, X3=X3, f1=f1, f2=f2, f3=f3, P1=P1, P2=P2, P3=P3)
    (I1+I2)
end


function prepare_containers_ltl(sources, ϵ, Nt, params::Constants, nmodel::N_Model)
    @unpack Δt, α, rb, kg, Δt̃ = params

    n = 10
    Ns = length(sources)
    # Evaluation points 
    # DO NOT USE FOR SELF-RESPONSE
    distances = zeros(Ns, Ns)
    NR = zeros(Int, Ns, Ns)
    K_min = zeros(Int, Ns, Ns)

    for j in 1:Ns
        for i in 1:j-1
            source = sources[i]
            target = sources[j]
            σ = FiniteLineSource.compute_distance_2D(source, target)
            setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H / 2 , σ=σ)
            N_r = compute_N_line(ϵ, setup, params)
            distances[i, j] = σ
            distances[j, i] = σ
            NR[i, j] = N_r
            NR[j, i] = N_r
        end
    end
    
    Nr = filter!(e -> e != 0, unique(NR))
    N = choose_blocks(Nr, Nt, p = 10)
    K = length(N) - 1

    for j in 1:Ns
        for i in 1:j-1
            Km = min(findlast(x -> x <= NR[i, j], N), K)
            K_min[i, j] = Km
            K_min[j, i] = Km
        end
    end
   
    ζ = zeros(0)
    W = zeros(0)
    indices = zeros(Int64, K+1)

    D_eff = sum([source.D for source in sources])/length(sources)
    H_eff = sum([source.H for source in sources])/length(sources)
    z_eff = D_eff + H_eff/2

    for i in 1:K
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ, n, D_eff, H_eff, z_eff, params, nmodel)
        @views indices[i+1:end] .+= N_block
    end

    @views ranges = [indices[i]+1:indices[i+1] for i in eachindex(indices[1:end-1])]
    @views Kranges = [index+1:indices[end] for index in indices[1:end-1]]

    #sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout = bakhalov_discretization(10, setup)
    
    # Preallocate objects
    F = zeros(length(ζ))
    expt = @. exp(-ζ^2*Δt̃)
    expNin = zeros(length(ζ))
    expNout = zeros(length(ζ))

    for i in 1:K
        @inbounds @. @views expNin[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i]*Δt̃)
        if i < K
            @inbounds @. @views expNout[ranges[i]] = exp(-ζ[ranges[i]]^2*N[i+1]*Δt̃)
        end
    end

    load_delays = [CircularBuffer{Float64}(N[i+1] - N[i]) for i in 1:K]
    load_buffer = CircularBuffer{Float64}(N[1])
    fill!(load_buffer, 0.)
    for load_delay in load_delays
        fill!(load_delay, 0.)
    end

    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    X, Fx, PP, XT = FiniteLineSource.create_bin_containers(bins)

    HM = zeros(length(ζ), length(sources), length(sources))

    C = 1 / (2π^2*kg)

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

        @inbounds HM[k, i, j] = C / target.H * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_double_line(lineparams, lineparamsT; x1=X[bin1], x2=X[bin2], x3=X[bin3], X1=XT[bin1], X2=XT[bin2], X3=XT[bin3], P1=PP[bin1], P2=PP[bin2], P3=PP[bin3], f1=Fx[bin1], f2=Fx[bin2], f3=Fx[bin3])
        @inbounds HM[k, j, i] = HM[k, i, j]
    end

    qin = zeros(length(ζ))
    qout = zeros(length(ζ))

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
        qin, 
        qout,
        #=sr_ζ,
        sr_w,
        sr_F,
        sr_expt,
        sr_expNout,
        sr_Ic,
        sr_Icout=#
    )
end