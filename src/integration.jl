function integrate(f, a, b, n)
    x, w = gausslegendre(n+1)
    m = (b-a)/2
    c = (b+a)/2
    @. x = x * m + c 
    m * dot(f.(x), w)
end

function integrate_x_w(f, a, b, n)
    x, w = gausslegendre(n+1)
    m = (b-a)/2
    c = (b+a)/2
    @. x = x * m + c
    @. w = m * w
    dot(f.(x), w), x, w
end


function adaptive_integration(f, a, b; tol = 10^-10, n = 3, min_interval = 10^-15)
    mid = (a+b)/2   

    I = integrate(f, a, b, n)
    I1, x1, w1 = integrate_x_w(f, a, mid, n)
    I2, x2, w2 = integrate_x_w(f, mid, b, n)

    I_sub = I1 + I2
    error = abs(I - I_sub)

    if error > tol && (b-a)/2 > min_interval
        In1, xn1, wn1 = adaptive_integration(f, a, mid, tol=tol, n=n, min_interval=min_interval)
        In2, xn2, wn2 = adaptive_integration(f, mid, b, tol=tol, n=n, min_interval=min_interval)
        return In1 + In2, [xn1; xn2], [wn1; wn2]
    else
        return I_sub, [x1; x2], [w1; w2]
    end
end
