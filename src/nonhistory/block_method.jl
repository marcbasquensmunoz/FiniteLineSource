
"""
Computes the blocks to be used
"""
function choose_blocks(Nr; p = 10)
    Nmin = minimum(Nr)
    Nmax = maximum(Nr)
    Ncurrent = Nmin
    N = zeros(Int, 0)

    while Ncurrent < 10*Nmax || length(N) <= 3
        push!(N, Int(floor(Ncurrent)))
        Ncurrent *= p
    end
    !(Nmin in N) && push!(N, Nmin)
    sort!(N)

    for (na, nb) in zip(N[1:end-1], N[2:end])
        if isempty(filter(n -> n in na:nb, Nr)) 
            deleteat!(N, findfirst(x->x==na, N))
        end
    end
    return N
end


"""
Compute the number of steps N that can be skipped for a line to point
"""
function compute_N_line(σ, D, H, z, ϵ, params::Constants) 
    @unpack Δt, α, kg = params
    f(N) = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z)^2) / sqrt(4α * Δt * N))/(4*π*kg*sqrt(σ^2 + (zp - z)^2)), D, D+H)[1] - ϵ
    Int(floor(find_zero(f, 10σ^2)))
end


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
function compute_ζ_points(N, No, ϵ, Q, n, params::Constants)
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
    sort!(segbuf, by=x->x.a)
    n_seg = length(segbuf)

    Nζ = n_seg*n

    ζ = zeros(Nζ)
    W = zeros(Nζ)
    x, w = gausslegendre(n)

    for (i, segment) in enumerate(segbuf)
        m = (segment.b-segment.a)/2
        c = (segment.b+segment.a)/2 
        @. ζ[(i-1)*n+1:i*n] = m*x + c
        @. W[(i-1)*n+1:i*n] = m * w
    end

    return ζ, W
end

function taylor(r, σ, m)
    f(x) = 1/(x^2-σ^2)^(1/2)
    # This returns the Taylor coefficients, to get the derivatives we need to multiply by m!
    derivatives(f, r, 1., Val(m)).partials
end

"""
Compute the minimum number of terms n to use in the Legendre expansion of the function 1/sqrt(x^2-σ^2)
to obtain an error lower than ϵ in the interval [a, b].
"""
function N_bound(ϵ, a, b, σ; C = 1, m = 30)
    #@show a, b
    M = 1:m
    Vm = C .* abs.(taylor(a, σ, m) .- taylor(b, σ, m))

    # Here we use the Stirling approximation for the factorial in order to avoid computing the factorials
    Int(ceil(minimum(@. (2π*M)^(1/(2M+1)) * (M/exp(1))^(M/(M + 1/2))*(Vm / (ϵ * sqrt(π*(M+0.5))))^(1/(M+0.5)) + M)))
end

compute_ζ_points_line(N, No, ϵ, Q, n, setup::SegmentToPoint, params::Constants) =
    compute_ζ_points_line(N, No, ϵ, Q, n, setup.D, setup.H, setup.z, params::Constants)

"""
Compute the nodes ζ and weights W suitable to integrate the function F after skipping N steps
"""
function compute_ζ_points_line(N, No, ϵ, Q, n, D, H, z, params::Constants)
    @unpack Δt, α, rb = params
    Δt̃ = Δt*α/rb^2

    heatwave(σ) = quadgk(zp -> erfc(sqrt(σ^2 + (zp - z)^2) / sqrt(4α * Δt * N))/(4*π*kg*sqrt(σ^2 + (zp - z)^2)), D, D+H)[1] - ϵ
    σ = find_zero(heatwave, sqrt(N))
    σ̃ = σ/rb

    a = 0.
    b = sqrt(-log(ϵ/Q) / (N*Δt̃))

    setup = SegmentToPoint(D=D, H=H, z=z_eval, σ=σ)

    Nl = Int(ceil(-log10(ϵ))) * 17 
    xb, wb = gausslegendre(Nl+1) 
    Pb = [wb[s] * Pl(xb[s], k) for k = 0:Nl, s = 1:Nl+1]
    Xb = zeros(Nl+1)
    aux = zeros(Nl+1)

    int_sin(ζ) = compute_kernel_line(ζ, rb, setup, x=xb, X=Xb, P=Pb, f=aux, atol=ϵ) 

    guide(ζ) = (exp(-ζ^2*N*Δt̃) + 50 * exp(-ζ^2*No*Δt̃)) * int_sin(ζ) * (1 - exp(-ζ^2*Δt̃)) / ζ
    _, _, segbuf = quadgk_segbuf(guide, a, b, order=n, atol=ϵ)
    sort!(segbuf, by=x->x.a)
    n_seg = length(segbuf)

    Nζ = n_seg*n

    ζ = zeros(Nζ)
    W = zeros(Nζ)
    x, w = gausslegendre(n)

    for (i, segment) in enumerate(segbuf)
        m = (segment.b-segment.a)/2
        c = (segment.b+segment.a)/2 
        @. ζ[(i-1)*n+1:i*n] = m*x + c
        @. W[(i-1)*n+1:i*n] = m * w
    end

    return ζ, W
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
    x1::Vector{T}
    x2::Vector{T}
    x3::Vector{T}
    X1::Vector{T}
    X2::Vector{T}
    X3::Vector{T}
    f1::Vector{T}
    f2::Vector{T}
    f3::Vector{T}
    P1::Matrix{T}
    P2::Matrix{T}
    P3::Matrix{T}
end

function LineKernelParams(setup::SegmentToPoint, ω, ϵ)
    @unpack D, H, z, σ = setup

    rB = sqrt(σ^2 + (z - D - H)^2 )
    rT = sqrt(σ^2 + (z - D)^2     )
    
    r1 = (D < z && z < D+H) ? σ : min(rB, rT)
    r2 = min(rB, rT)
    r3 = max(rB, rT)
    #rs = r1 + 4π/ω
    rs = r1 + min(4π/ω, (r2-r1)) 
    n1, n2 = 0, 0

    #@show r1, rs, r2, r3
    if r1 != r2
        #@info "Computing n corresponding to I2"
        n1 = N_bound(ϵ / (6 * sqrt(r2-rs)), rs, r2, σ, C=2) + 20
        #@show rs, r2, σ, ϵ / (6 * sqrt(r2-rs)), n1
    end

    if r2 != r3
        #@info "Computing n corresponding to I1"
        rl = max(r2, rs)
        n2 = N_bound(ϵ / (3 * sqrt(r3-rl)), rl, r3, σ) + 20
        #@show rl, r3, σ, ϵ / (3 * sqrt(r3-rl)), n2
    end

    LineKernelParams(r1=r1, r2=r2, r3=r3, rs=rs, ω=ω, σ=σ, n1=n1, n2=n2, ϵ=ϵ)
end


function LineToLineKernelParams(setup::SegmentToSegment, ω, ϵ)
    @unpack D1, H1, D2, H2, σ = setup

    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     ) 
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     ) 
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : σ
    r2 = H2 > H1 ? (D2 > D1 ? rLL : σ) : (D2 + H2 > D1 + H1 ? rUR : σ)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : σ) : (D2 > D1 ? rLL : σ)
    r4 = D2 + H2 > D1 ? rUL : σ

    rs = r1 + 4π/ω

    α1 = D1 - D2 + H1
    α2 = min(H1, H2)
    α3 = D2 - D1 + H2

    x1, x2, x3 = zeros(0), zeros(0), zeros(0)
    X1, X2, X3 = zeros(0), zeros(0), zeros(0)
    f1, f2, f3 = zeros(0), zeros(0), zeros(0)
    P1, P2, P3 = zeros(0, 0), zeros(0, 0), zeros(0, 0)

    if r1 != r2
        n1 = N_bound(ϵ / (3α1 * sqrt(r2-rs)), rs, r2, σ)
        x1, w1 = gausslegendre(n1+1) 
        P1 = [w1[s] * Pl(x1[s], k) for k in 0:n1, s in 1:n1+1]
        X1 = zeros(n1+1)
        f1 = zeros(n1+1)
    end

    if r2 != r3
        rl = max(r2, rs)
        n2 = N_bound(ϵ / (3α2 * sqrt(r3-rl)), rl, r3, σ)
        x2, w2 = gausslegendre(n2+1) 
        P2 = [w2[s] * Pl(x2[s], k) for k in 0:n2, s in 1:n2+1]
        X2, f2 = zeros(n2+1), zeros(n2+1)
    end

    if r3 != r4
        rl = max(r3, rs)
        n3 = N_bound(ϵ / (3α3 * sqrt(r4-rl)), rl, r4, σ)
        x3, w3 = gausslegendre(n3+1) 
        P3 = [w3[s] * Pl(x3[s], k) for k in 0:n3, s in 1:n3+1]
        X3, f3 = zeros(n3+1), zeros(n3+1)
    end

    LineToLineKernelParams(r1=r1, r2=r2, r3=r3, r4=r4, rs=rs, α1=α1, α2=α2, α3=α3, ω=ω, σ=σ, x1=x1, x2=x2, x3=x3, X1=X1, X2=X2, X3=X3, f1=f1, f2=f2, f3=f3, P1=P1, P2=P2, P3=P3)
end

function compute_kernel_line(params::LineKernelParams; x1=nothing, x2=nothing, X1=nothing, X2=nothing, P1=nothing, P2=nothing, f1=nothing, f2=nothing)
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
        m1 = (r2-rs)/2
        c1 = (r2+rs)/2
        @. X1 = m1*x1 + c1

        @. f1 = 2 / sqrt(X1^2-σ^2)
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)

        # This version creates less allocations, but runs slower
        #I_osc = sqrt(m*π/(2ω)) * dot(X, P, f)
        I_osc += sqrt(m1*π/(2ω)) * X1' * P1 * f1
    end

    if r2 != r3    
        #@info "Computing I1"
        rl = r1 == σ ? max(rs, r2) : r2
        m2 = (r3-rl)/2
        c2 = (r3+rl)/2
        @. X2 = m2*x2 + c2

        @. f2 = 1. / sqrt(X2^2-σ^2)
        besselj!(X2, 1/2:(n2+1/2), m2*ω)
        @. X2 = X2 * imag(exp(im*ω*c2) * im^(0:n2)) * (2(0:n2)+1)
        I_osc += sqrt(m2*π/(2ω)) * X2' * P2 * f2
    end

    I_div + I_osc
end


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


function compute_kernel_double_line_part(params::LineToLineKernelParams; atol=1e-8)
    @unpack r1, r2, r3, r4, rs, ω, α1, α2, α3, σ, x1, x2, x3, X1, X2, X3, f1, f2, f3, P1, P2, P3 = params

    h(r) = r < r2 ? α1 : (r < r3 ? α2 : α3)
    I_div = 0.
    I_osc = 0.

    if r1 == σ
        mult = r1 == r2 ? (r2 == r3 ? α3 : α2) : α1
        C = mult * sin(σ*ω)
        h_reg(r) = r == σ ? 0. : (sin(r*ω) * h(r) - C) / sqrt(r^2-σ^2)
        I1, E, count = quadgk_count(h_reg, r1, rs, atol=atol/3)
        I2 = C * acoth(rs/sqrt(rs^2-r1^2))
        I_div = I1 + I2
        #@show I_div
    end

    if r1 != r2
        m1 = (r2-rs)/2
        c1 = (r2+rs)/2
        @. X1 = m1*x1 + c1
        n1 = length(X1)-1

        @. f1 = α1 / sqrt(X1^2-σ^2)
        besselj!(X1, 1/2:(n1+1/2), m1*ω)
        @. X1 = X1 * imag(exp(im*ω*c1) * im^(0:n1)) * (2(0:n1)+1)
        #@show sqrt(m1*π/(2ω)) * X1' * P1 * f1
        I_osc += sqrt(m1*π/(2ω)) * X1' * P1 * f1
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
        #@show sqrt(m2*π/(2ω)) * X2' * P2 * f2
        I_osc += sqrt(m2*π/(2ω)) * X2' * P2 * f2
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
        #@show sqrt(m3*π/(2ω)) * X3' * P3 * f3
        I_osc += sqrt(m3*π/(2ω)) * X3' * P3 * f3
    end
    I_osc_linear = (cos(ω*r1) - cos(ω*r2) - cos(ω*r3) + cos(ω*r4))/ω
    #@show I_osc_linear
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
    @. f = h1(X) / real(sqrt(complex(X)^2-σ^2))

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
