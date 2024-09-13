using DataStructures
function integrate(f, a, b, n)
    x, w = gausslegendre(n+1)
    m = (b-a)/2
    c = (b+a)/2
    @. x = x * m + c 
    m * dot(f.(x), w)
end

@with_kw struct IntegrationSegment{T <: Number} @deftype T
    a
    b
    I
    E
end
Base.isless(x::IntegrationSegment{T}, y::IntegrationSegment{T}) where T <: Number = x.E < y.E

function adaptive_gk_segments(f, a::T, b::T; rtol = sqrt(eps())) where T <: Number
    n_pre = Int(floor((b-a) / 2))
    n = max(n_pre%2 == 0 ? n_pre : n_pre+1, 8)
    x, w, gw = QuadGK.kronrod(n)
    heap = MutableBinaryMaxHeap{IntegrationSegment{T}}()

    I, E = evalrule(f, a, b, x, w, gw)
    push!(heap, IntegrationSegment(a, b, I, E))

    while E > rtol * abs(I)
        s = pop!(heap)
        mid = (s.a+s.b) * convert(eltype(a), 0.5)
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I =+ I1 + I2 - s.I
        E =+ E1 + E2 - s.E
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
    end

    extract_all!(heap)
end

function adaptive_gk(f, a::T, b::T; rtol = sqrt(eps()), atol = eps()) where T <: Number
    x, w, gw = QuadGK.kronrod(8)
    heap = MutableBinaryMaxHeap{IntegrationSegment{T}}()

    I, E = evalrule(f, a, b, x, w, gw)
    push!(heap, IntegrationSegment(a, b, I, E))

    Ns = 1

    while E > rtol * abs(I) && abs(E) > atol
        s = pop!(heap)
        mid = (s.a+s.b) * convert(eltype(a), 0.5)
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I =+ I1 + I2 - s.I 
        E =+ E1 + E2 - s.E 
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
        Ns += 1
    end
    extract_nodes_weights(heap, x, w, Ns)
end

function extract_nodes_weights(heap::MutableBinaryMaxHeap{IntegrationSegment{T}}, x, w, Ns) where T <: Number
    n = 2*length(x)-1
    res_x = zeros(T, Ns*n)
    res_w = zeros(T, Ns*n)

    i = 1

    while !isempty(heap)
        s = pop!(heap)
        @views rescale!(s.a, s.b, x, w, res_x[(i-1)*n+1:i*n], res_w[(i-1)*n+1:i*n])

        i += 1
    end

    res_x, res_w
end

function rescale!(a, b, x, w, X, W)
    n = length(x)
    m = (b-a) * 0.5
    c = (b+a) * 0.5
    @. X[1:n] = m * x + c
    @views @. X[n+1:end] = -m * x[1:end-1] + c
    @. W[1:n] = m * w
    @views @. W[n+1:end] = m * w[1:end-1]
end

function evalrule(f, a::T, b::T, x::Vector{T}, w::Vector{S}, gw::Vector{S}) where {T <: Number, S <: Number}
    s = (b-a) * convert(eltype(a), 0.5)
    
    Fg = @views @. f(a + (1+x[2:2:end])*s) + f(a + (1-x[2:2:end])*s)
    Fk = @views @. f(a + (1+x[1:2:end-1])*s) + f(a + (1-x[1:2:end-1])*s)

    Ig = dot(Fg, gw)
    @views Ik = f(a + s) * w[end] + dot(Fg, w[2:2:end]) + dot(Fk, w[1:2:end-1])

    return Ik * s, abs((Ik - Ig) * s)
end
