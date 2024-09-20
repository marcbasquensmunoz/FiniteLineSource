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

function adaptive_gk_segments(f, a::T, b::T; rtol = sqrt(eps()), atol = eps()) where T <: Number
    n_pre = Int(floor((b-a) / 2))
    n = max(n_pre%2 == 0 ? n_pre : n_pre+1, 8)
    x, w, gw = QuadGK.kronrod(n)
    heap = MutableBinaryMaxHeap{IntegrationSegment{T}}()

    eval1, eval2 = zeros(length(gw)), zeros(length(gw))

    A = @MVector zeros(2)
    B = @MVector zeros(2)
    C = @MVector zeros(2)

    evalrule!(A, f, a, b, x, w, gw, eval1, eval2)
    push!(heap, IntegrationSegment(a, b, A[1], A[2]))

    while A[2] > rtol * abs(A[1]) && abs(A[2]) > atol
        s = pop!(heap)
        mid = (s.a+s.b) * convert(eltype(a), 0.5)
        evalrule!(B, f, s.a, mid, x, w, gw, eval1, eval2)
        evalrule!(C, f, mid, s.b, x, w, gw, eval1, eval2)
        A[1] += B[1] + C[1] - s.I 
        A[2] += B[2] + C[2] - s.E 
        push!(heap, IntegrationSegment(s.a, mid, B[1], B[2]))
        push!(heap, IntegrationSegment(mid, s.b, C[1], C[2]))
    end

    extract_all!(heap)
end

function adaptive_gk(f, a, b; rule = QuadGK.kronrod(8), rtol = sqrt(eps()), atol = eps()) 
    x, w, gw = rule
    heap = MutableBinaryMaxHeap{IntegrationSegment{typeof(a)}}()
        #MutableBinaryHeap{IntegrationSegment{T}, DataStructures.FasterReverse}()

    eval1, eval2 = zeros(length(gw)), zeros(length(gw))

    A = @MVector zeros(2)
    B = @MVector zeros(2)
    C = @MVector zeros(2)

    evalrule!(A, f, a, b, x, w, gw, eval1, eval2)
    push!(heap, IntegrationSegment(a, b, A[1], A[2]))

    Ns = 1

    while A[2] > rtol * abs(A[1]) && abs(A[2]) > atol
        s = pop!(heap)
        mid = (s.a+s.b) * convert(typeof(s.a), 0.5)

        evalrule!(B, f, s.a, mid, x, w, gw, eval1, eval2)
        evalrule!(C, f, mid, s.b, x, w, gw, eval1, eval2)
        A[1] += B[1] + C[1] - s.I 
        A[2] += B[2] + C[2] - s.E 
        push!(heap, IntegrationSegment(s.a, mid, B[1], B[2]))
        push!(heap, IntegrationSegment(mid, s.b, C[1], C[2]))
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

function evalrule!(I, f, a, b, x, w, gw, aux1, aux2)
    s = (b-a) * convert(typeof(a), 0.5)

    @. @views aux1 = f(a + (1+x[2:2:end])*s) + f(a + (1-x[2:2:end])*s)
    @. @views aux2 = f(a + (1+x[1:2:end-1])*s) + f(a + (1-x[1:2:end-1])*s)

    Ig = dot(aux1, gw)
    @views Ik = f(a + s) * w[end] + dot(aux1, w[2:2:end]) + dot(aux2, w[1:2:end-1])

    I[1] = Ik * s
    I[2] = abs((Ik - Ig) * s)

    return nothing
end
