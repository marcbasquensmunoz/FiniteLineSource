mutable struct Buffer{T <: Number}
    segments::Vector{QuadGK.Segment{T, T, T}}
    X::Vector{T}
    W::Vector{T}
end
Buffer(segments) = Buffer(segments, zeros(0), zeros(0))

function rescale!(a, b, x, w, X, W)
    n = length(x)
    m = (b-a) * 0.5
    c = (b+a) * 0.5
    @. X[1:n] = m * x + c
    @views @. X[n+1:end] = -m * x[1:end-1] + c
    @. W[1:n] = m * w
    @views @. W[n+1:end] = m * w[1:end-1]
    return nothing
end
function adaptive_nodes_and_weights(f, a, b; n = 8, buffer = nothing, rtol = sqrt(eps()))
    m = 2*n+1
    if isnothing(buffer)
        _, _, segs = quadgk_segbuf(f, a, b, order = n, rtol = rtol)
        buffer = Buffer(segs, zeros(m*length(segs)), zeros(m*length(segs)))
    else
        quadgk(f, a, b, order = n, rtol = rtol, segbuf = buffer.segments, eval_segbuf = buffer.segments)
        if length(buffer.segments) * m != length(buffer.X)
            buffer.X = zeros(m*length(buffer.segments))
            buffer.W = zeros(m*length(buffer.segments))
        end
    end
    x, w, _ = QuadGK.cachedrule(Float64, n)
    for (i, segment) in enumerate(buffer.segments)
        @views rescale!(segment.a, segment.b, x, w, buffer.X[(i-1)*m+1:i*m], buffer.W[(i-1)*m+1:i*m])
    end
    return buffer.X, buffer.W
end