
struct BlockMethod{T <: Number}
    ζ::Vector{T}
    F::Vector{T}
    expt::Vector{T}
    expNin::Vector{T}
    expNout::Vector{T}
    HM::Array{T, 3}
    load_delays::Vector{CircularBuffer{T}}
    load_buffer::CircularBuffer{T}
    ranges::Vector{UnitRange{Int}}
    Kranges::Vector{UnitRange{Int}}
    K_min::Matrix{Int}
    qinaux::Vector{T}
    qoutaux::Vector{T}
end

"""
Computes the blocks to be used
"""
function choose_blocks(Nr, Nt; p = 10)
    Nmin = minimum(Nr)
    Ncurrent = Nmin
    N = zeros(Int, 0)

    while Ncurrent < Nt
        push!(N, Int(floor(Ncurrent)))
        Ncurrent *= p
    end
    !(Nmin in N) && push!(N, Nmin)
    sort!(N)
    push!(N, Nt)

    @views for (na, nb) in zip(N[1:end-1], N[2:end])
        if isempty(filter(n -> n in na:nb, Nr)) 
            deleteat!(N, findfirst(x->x==na, N))
        end
    end
    return N
end


function taylor(r, σ, m)
    f(x) = 1/sqrt(x^2-σ^2)
    # This returns the Taylor coefficients, to get the derivatives we need to multiply by m!
    derivatives(f, r, 1., Val(m)).partials
end

"""
Compute the minimum number of terms n to use in the Legendre expansion of the function 1/sqrt(x^2-σ^2)
to obtain an error lower than ϵ in the interval [a, b].
"""
function N_bound(ϵ, a, b, σ, nmodel; C = 1., m = 30, n_max = 500)
    return Int(ceil(eval(nmodel, ϵ, σ, a, b)))
    M = 1:m
    #Vm = C .* abs.(taylor(a, σ, m) .- taylor(b, σ, m))
    Vm = abs.(taylor(a, σ, m) .- taylor(b, σ, m)) .* (M/exp(1)) .^ M .* sqrt.(2*π .* M) .* ((b-a)/2) .^ M
    # Here we use the Stirling approximation for the factorial in order to avoid computing the factorials
    n = Int(ceil(minimum(@. (2π*M)^(1/(2M+1)) * (Vm / (ϵ * sqrt(π*(M+0.5))))^(1/(M+0.5)) + M)))
    min(n_max, n)
end

function N_bound_1(ϵ, a, b, σ)
    #C = 0.3
    C = 0.65
    #α = 0.3386
    α = 0.116
    return Int(ceil(-log(ϵ/C)/α)) 
    Int(ceil(- 0.8 * log10(ϵ/5) * sqrt(b-a)/(a-σ)^(1/2)))
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

@with_kw struct LineKernelContainers{T <: Number} @deftype T
    x::Vector{T}
    X::Vector{T}
    f::Vector{T}
    P::Matrix{T}
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
    ϵ
    n1::Int
    n2::Int
    n3::Int
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
        n1 = N_bound(ϵ / (6), rs, r2, σ, nmodel, n_max=n_max)
        #@show rs, r2, σ, ϵ / (6 * sqrt(r2-rs)), n1
    end

    if r2 != r3
        #@info "Computing n corresponding to I1"
        rl = max(r2, rs)
        n2 = N_bound(ϵ / (3 * sqrt(r3-rl)), rl, r3, σ, nmodel, n_max=n_max)
        #@show rl, r3, σ, ϵ / (3 * sqrt(r3-rl)), n2
    end

    LineKernelParams(r1=r1, r2=r2, r3=r3, rs=rs, ω=ω, σ=σ, n1=n1, n2=n2, ϵ=ϵ)
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

function create_bin_containers(bins)
    X = Vector{Vector{Float64}}()
    PP = Vector{Matrix{Float64}}()
    F = Vector{Vector{Float64}}()
    XT = Vector{Vector{Float64}}()

    for n in bins
        x, w = gausslegendre(n+1)
        P = zeros(n+1, n+1)
        for s in 1:n+1
            @inbounds @views collectPl!(P[:, s], x[s], lmax=n)
            @inbounds @views @. P[:, s] *= w[s]
        end
        f = zeros(n+1)
        xt = zeros(n+1)

        push!(X, x)
        push!(F, f)
        push!(PP, P)
        push!(XT, xt)
    end

    X, F, PP, XT
end
get_bin(n, bins) = findfirst(m -> n < m, bins)

function evolve!(I, q, block::BlockMethod{T}) where {T <: Number}
    @unpack ζ, F, expt, expNin, expNout, HM, load_delays, load_buffer, 
        ranges, Kranges, K_min, qinaux, qoutaux#=, sr_ζ, sr_w, sr_F, sr_expt, sr_expNout, sr_Ic, sr_Icout =#= block

    for (nt, qt) in enumerate(q)
        current_q = load_buffer[1]
        push!(load_buffer, qt)
        
        #=
        @. sr_F = sr_expt * (sr_F - qt / sr_ζ + sr_expNout * current_q / sr_ζ)
        for target in 1:size(block.K_min)[1]
            I[target, nt] += dot(sr_F, sr_w) + qt * sr_Ic - current_q * sr_Icout
        end
        @. sr_F = sr_F + (qt - sr_expNout * current_q) / sr_ζ
        =#
        for i in eachindex(load_delays)
            qin = current_q
            qout = i == length(load_delays) ? 0. : load_delays[i][1]
            push!(load_delays[i], qin)
            current_q = qout
            @. qinaux[ranges[i]] = qin
            @. qoutaux[ranges[i]] = qout
        end

        @. F = expt * F + (qinaux * expNin - qoutaux * expNout)

        bh_indices = 1:size(block.K_min)[1]
        for target in bh_indices
            for source in bh_indices
                if source == target continue end
                range = Kranges[K_min[source, target]]
                @views I[target, nt] += dot(F[range], HM[range, target, source])
            end
        end
    end
end
