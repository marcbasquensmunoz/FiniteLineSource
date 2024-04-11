function integrate(f, a, b, n)
    x, w = gausslegendre(n+1)
    m = (b-a)/2
    c = (b+a)/2
    @. x = x * m + c
    m * dot(f.(x), w)
end
