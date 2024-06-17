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
        mid = (s.a+s.b) * convert(eltype(x), 0.5)
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I = I - s.I + I1 + I2
        E = E - s.E + E1 + E2
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
    end

    extract_all!(heap)
end

function adaptive_gk(f, a, b; rtol = sqrt(eps()))
    x, w, gw = QuadGK.kronrod(8)
    heap = MutableBinaryMaxHeap{IntegrationSegment{Float64}}()

    I, E = evalrule(f, a, b, x, w, gw)
    push!(heap, IntegrationSegment(a, b, I, E))

    while E > rtol * abs(I)
        s = pop!(heap)
        mid = (s.a+s.b) * 0.5
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I = I - s.I + I1 + I2
        E = E - s.E + E1 + E2
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
    end
    extract_nodes_weights(heap, x, w)
end

function extract_nodes_weights(heap, x, w)
    res_x = zeros(0)
    res_w = zeros(0)

    while !isempty(heap)
        s = pop!(heap)
        xs, ws = rescale_x_w(s.a, s.b, x, w)
        append!(res_x, xs)
        append!(res_w, ws)
    end

    res_x, res_w
end

function rescale_x_w(a, b, x, w)
    m = (b-a) * 0.5
    c = (b+a) * 0.5
    m .* [x; -x[1:end-1]] .+ c, m .* [w; w[1:end-1]]
end

function evalrule(f, a, b, x, w, gw)
    s = (b-a) * 0.5

    Ik = f(a + s) * w[end]
    Ig = zero(Ik)

    for i in eachindex(gw)
        fg = f(a + (1+x[2i])*s) + f(a + (1-x[2i])*s)
        fk = f(a + (1+x[2i-1])*s) + f(a + (1-x[2i-1])*s)
        Ig += fg * gw[i]
        Ik += fg * w[2i] + fk * w[2i-1]
    end
    Ik_s, Ig_s = Ik * s, Ig * s
    E = abs(Ik_s - Ig_s)
    return Ik_s, E
end
