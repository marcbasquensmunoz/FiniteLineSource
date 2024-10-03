
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
    if isnothing(buffer)
        _, _, segs = quadgk_segbuf(f, a, b, order = n, rtol = rtol)
        buffer = segs
    else
        quadgk(f, a, b, order = n, rtol = rtol, segbuf = buffer, eval_segbuf = buffer)
    end
    x, w, _ = QuadGK.cachedrule(Float64, n)
    m = 2*n+1
    X, W = zeros(m*length(buffer)), zeros(m*length(buffer))
    for (i, segment) in enumerate(buffer)
        @views rescale!(segment.a, segment.b, x, w, X[(i-1)*m+1:i*m], W[(i-1)*m+1:i*m])
    end
    return X, W
end