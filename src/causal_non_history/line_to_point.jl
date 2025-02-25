
compute_distance_2D(s, t) = sqrt((s.x - t.x)^2 + (s.y - t.y)^2)

struct N_Model{M, T <: Number}
    model::M
    mins::Vector{T}
    norms::Vector{T}
    aux::Vector{Float32}
end
function eval(nmodel, ϵ, σ, a, b)
    @unpack model, mins, norms, aux = nmodel
    aux[1] = (-log10(Float32(ϵ)) - mins[1]) / norms[1]
    aux[2] = (log(Float32(σ)) - mins[2]) / norms[2]
    aux[3] = (log(Float32(a)) - mins[3]) / norms[3]
    aux[4] = (log(Float32(b)) - mins[4]) / norms[4]
    Float64(model(aux)[1])
end

function load_nmodel()
    k = 32
    model = Chain(
        Dense(4 => k, relu),   
        Dense(k => k, relu),   
        Dense(k => 1)
    )

    mins = Float32[4.0, 0.0, 0.6931471824645996, 0.6931471824645996]
    norms = Float32[8.0, 6.907755374908447, 10.819777965545654, 13.122363567352295] 
    model_state = JLD2.load("N_bound.jld2", "model_state")
    Flux.loadmodel!(model, model_state)
    N_Model(model, mins, norms, zeros(Float32, 4))
end

@with_kw struct LineKernelParams{T <: Number} @deftype T
    r1
    r2
    r3
    rs
    ω
    σ
    ϵ
    n1::Int
    n2::Int
end

function LineKernelParams(setup::SegmentToPoint, ω, ϵ, nmodel; n_max=500)
    @unpack D, H, z, σ = setup

    rB = sqrt(σ^2 + (z - D - H)^2 )
    rT = sqrt(σ^2 + (z - D)^2     )
    
    r1 = (D < z && z < D+H) ? σ : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)
    #rs = r1 + 4π/ω
    #rs = r1 + min(4π/ω, (r2-r1)*0.1) 
    rs = min(r1 + max(2, (r2-r1) * 0.01), r2)
    n1, n2 = 0, 0

    if r1 != r2
        #@info "Computing n corresponding to I2"
        n1 = N_bound(ϵ / (6), rs, r2, σ, nmodel)
        #@show rs, r2, σ, ϵ / (6 * sqrt(r2-rs)), n1
    end

    if r2 != r3
        #@info "Computing n corresponding to I1"
        rl = max(r2, rs)
        n2 = N_bound(ϵ / (3 * sqrt(r3-rl)), rl, r3, σ, nmodel)
        #@show rl, r3, σ, ϵ / (3 * sqrt(r3-rl)), n2
    end

    LineKernelParams(r1=r1, r2=r2, r3=r3, rs=rs, ω=ω, σ=σ, n1=n1, n2=n2, ϵ=ϵ)
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
function compute_ζ_points_line!(ζ, W, N, No, ϵ, n, D, H, z, params::Constants, nmodel)
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

    Nl = Int(ceil(-log10(ϵ))) * 20
    xb, wb = gausslegendre(Nl+1) 
    Pb = zeros(Nl+1, Nl+1)
    for s in 1:Nl+1
        @inbounds @views collectPl!(Pb[:, s], xb[s], lmax=Nl)
        @inbounds @views @. Pb[:, s] *= wb[s]
    end
    BV = LineKernelContainers(x=xb, P=Pb)

    int_sin(ζ) = compute_kernel_line(LineKernelParams(setup, ζ/rb, ϵ, nmodel), BV1=BV, BV2=BV)
    #int_sin(ζ) = compute_kernel_line(ζ, rb, setup, x=xb, X=Xb, P=Pb, f=aux, atol=ϵ) 

    guide(ζ) = (No-N) * (exp(-ζ^2*N*Δt̃) + exp(-ζ^2*No*Δt̃)) * int_sin(ζ) * (1 - exp(-ζ^2*Δt̃)) / ζ
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

function compute_kernel_line(params::LineKernelParams; BV1, BV2)#x1=nothing, x2=nothing, X1=nothing, X2=nothing, P1=nothing, P2=nothing, f1=nothing, f2=nothing)
    @unpack r1, r2, r3, rs, ω, σ, ϵ, n1, n2 = params

    h(r) = r < r2 ? 2. : 1.
    I_div = 0.
    I_osc = 0.

    if r1 == σ
        #@info "Computing Is"
        mult = r1 == r2 ? 1. : 2.
        C = mult * sin(σ*ω)
        h_reg(r) = r == σ ? 0. : (sin(r*ω) * h(r) - C) / sqrt(r^2-σ^2)
        # We could do the integration vectorized in ω
        I1, _ = quadgk(h_reg, r1, rs, atol=ϵ/3)
        I2 = C * acoth(rs/sqrt(rs^2-r1^2))
        I_div = I1 + I2
    end

    if r1 != r2
        #@info "Computing I2"
        x1 = BV1.x
        X1 = BV1.X
        P1 = BV1.P
        f1 = BV1.f
        n1 = BV1.n

        m1 = (r2-rs)/2
        c1 = (r2+rs)/2
        @. X1 = m1*x1 + c1

        @. f1 = 2. / sqrt(abs(X1^2-σ^2))
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        I_osc = sqrt(m1*π/(2ω)) * dot(X1, P1, f1)

        #=
        rs1 = rs + 1.

        m1 = (rs1-rs)/2
        c1 = (rs1+rs)/2
        @. X1 = m1*x1 + c1
        @. f1 = 2. / sqrt(abs(X1^2-σ^2))
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        I_osc += sqrt(m1*π/(2ω)) * X1' * P1 * f1

        m1 = (r2-rs1)/2
        c1 = (r2+rs1)/2
        @. X1 = m1*x1 + c1
        @. f1 = 2. / sqrt(abs(X1^2-σ^2))
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        I_osc += sqrt(m1*π/(2ω)) * X1' * P1 * f1
        =#
    end

    if r2 != r3    
        x2 = BV2.x
        X2 = BV2.X
        P2 = BV2.P
        f2 = BV2.f
        n2 = BV2.n

        #@info "Computing I1"
        rl = r1 == σ ? max(rs, r2) : r2
        m2 = (r3-rl)/2
        c2 = (r3+rl)/2
        @. X2 = m2*x2 + c2

        @. f2 = 1. / sqrt(X2^2-σ^2)
        besselj!(X2, 1/2:(n2+1/2), m2*ω)
        @. X2 = X2 * imag(exp(im*ω*c2) * im^(0:n2)) * (2(0:n2)+1)
        I_osc += sqrt(m2*π/(2ω)) * dot(X2, P2, f2)
    end

    I_div + I_osc
end

#=
function compute_kernel_line(ζ, rb, setup::SegmentToPoint; x, X, P, f, atol=1e-8)
    @unpack D, H, z, σ = setup
    rmin = σ
    rB = sqrt(σ^2 + (z - D - H)^2 )
    rT = sqrt(σ^2 + (z - D)^2     )
    
    r1 = (D < z && z < D+H) ? rmin : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)

    ω = ζ/rb
    h(r) = r < r2 ? 2. : 1.

    split = r1
    I_div = 0.

    if r1 == σ
        split = r1 + min(4π/ω, (r3-r1)*0.1)
        mult = r1 == r2 ? 1. : 2.
        C = mult * sin(σ*ω)
        h_reg(r) = r == σ ? 0. : (sin(r*ω) * h(r) - C) / sqrt(r^2-σ^2)
        # We could do the integration vectorized in ω
        I1, E, count = quadgk_count(h_reg, r1, split, atol=atol)
        #@show E, count
        I2 = C * acoth(split/sqrt(split^2-r1^2))
        I_div = I1 + I2
    end
    #@show split, r3

    m = (r3-split)/2
    c = (r3+split)/2

    n = length(x)-1
    @. X = m*x + c
    @. f = h(X) / sqrt(X^2-σ^2)

    besselj!(X, 1/2:(n+1/2), m*ω)
    @. X = X * imag(exp(im*ω*c) * im^(0:n)) * (2(0:n)+1)

    # This version creates less allocations, but runs 6 μs slower
    #I_osc = sqrt(m*π/(2ω)) * dot(X, P, f)
    I_osc = sqrt(m*π/(2ω)) * X' * P * f
    I_div + I_osc
end
=#


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

function compute_distance(::SegmentToPoint, source, target, params, ϵ)
    σ = compute_distance_2D(source, target)
    setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H / 2 , σ=σ)
    N_r = compute_N_line(ϵ, setup, params)
    σ, N_r
end

function compute_ζ_discretization!(ζ, W, indices, ::SegmentToPoint; sources, N, n, ϵ, constants, nmodel)
    K = length(N) - 1
    D_eff = sum([source.D for source in sources])/length(sources)
    H_eff = sum([source.H for source in sources])/length(sources)
    z_eff = D_eff + H_eff/2
    for i in 1:K
        N_block = compute_ζ_points_line!(ζ, W, N[i], N[i+1], ϵ, n, D_eff, H_eff, z_eff, constants, nmodel)
        @views indices[i+1:end] .+= N_block
    end
end

function compute_H!(HM, ::SegmentToPoint; ζ, W, expt, sources, distances, constants::Constants, ϵ, nmodel)
    @unpack kg, rb = constants
    C = 1 / (2π^2*kg)

    bins = [10, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250]
    containers = create_bin_containers(bins)

    for j in eachindex(sources), i in 1:j-1, (k, ζζ) in enumerate(ζ)
        σ = i == j ? rb : distances[i, j]
        source = sources[j]
        target = sources[i]
        setup = SegmentToPoint(D=source.D, H=source.H, z=target.D + target.H/2, σ=σ)
        lineparams = LineKernelParams(setup, ζζ/rb, ϵ, nmodel)
        bin1 = FiniteLineSource.get_bin(lineparams.n1, bins)
        bin2 = FiniteLineSource.get_bin(lineparams.n2, bins)
        @inbounds HM[k, i, j] = C * W[k] * (1 - expt[k]) / ζ[k] * compute_kernel_line(lineparams; BV1=containers[bin1], BV2=containers[bin2])
        @inbounds HM[k, j, i] = HM[k, i, j]
    end
end
