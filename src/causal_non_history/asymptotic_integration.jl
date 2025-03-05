@with_kw struct AsymptoticContainers{T<:Number}
    A::Matrix{T}
    v1::Vector{T}
    v2::Vector{T}
    w1::Vector{T}
    w2::Vector{T}
    aux::Vector{T}
    N::Int
end
function AsymptoticContainers(N::Int)
    if N%2 == 0 N += 1 end
    A = zeros(N+1, N+1)
    A[1, 1] = 1.
    
    for i in 2:N+1
        if i%2 == 0
            k = div(i, 2)
            p = k-1:2*(k-1)
            @. @views A[i, p+1] = - A[i-1, p+1] * (2p+1)
        else
            k = div(i-1, 2)
            A[i, k+1] = A[i-2, k] * (2k-1) * 2k
            A[i, 2k+1] = A[i-2, 2k-1] * (4k-3) * (4k-1)
            p = k+1:2k-1
            @. @views A[i, p+1] = A[i-2, p-1] * (2p-3)*(2p-1) + A[i-2, p] * (2p-1) * 2p
        end
    end
    AsymptoticContainers(A=A, v1=zeros(N+1), v2=zeros(N+1), w1=zeros(N+1), w2=zeros(N+1), aux=zeros(N+1), N=N)
end

function asymptotic_F!(x, containers::AsymptoticContainers)
    @unpack A, v1, v2, aux, N = containers
    r = 1/(x^2-1)
    sqrtr = sqrt(r)
    rp = 1.
    for p in 0:N-1
        #=
        v1[p+1] = r^(p+1/2) 
        v2[p+1] = p%2 == 1 ? x*r : 1.
        =#
        v1[p+1] = rp*sqrtr
        v2[p+1] = p%2 == 1 ? x*r : 1.
        rp *= r
    end
    mul!(aux, A, v1)
    aux .*= v2
    return nothing
end

function asymptotic_error(a, b, σ, ω, containers::AsymptoticContainers)
    @unpack N, A = containers

    σNI = let A=A, N=N
        x -> begin 
            res = 0.
            r = 1/(x^2-1)
            rp = 1.
            sqrtr = sqrt(r)
            for p in 0:N
                #res -= A[N+1, p+1] * r^(p+1) / (2p+1) 
                res -= A[N+1, p+1] * rp*sqrtr / (2p+1)
                rp *= r
            end
            return res 
        end
    end

    return abs(σNI(b/σ) - σNI(a/σ)) / (ω*σ)^N 

    #=

    σN = let A=A, N=N
        x -> begin 
        res = 0.
        for p in 0:N
            res += A[N+1, p+1] / sinh(x)^(2p) 
        end
        res *= N%2 == 1 ? cosh(x)/sinh(x)^2 : 1.
        return res 
    end
    end
    # Second best
    #return abs(quadgk(σN, acosh(a/σ), acosh(b/σ))[1]) / (ω*σ)^N 

    # Accurate
    f = let N=N, ω=ω, σ=σ
        x -> σN(x) * imag(exp(im*ω*σ*cosh(x))/(-im)^N)
    end

    1 / (ω*σ)^N * quadgk(f, acosh(a/σ), acosh(b/σ))[1]
    =#
end

function asymptotic(a, b, σ, ω, containers::AsymptoticContainers)
    @unpack N, w1, w2, aux = containers
    expa = exp(im*ω*a)
    expb = exp(im*ω*b)
    r = im/(ω*σ)
    rp = r
    for i in 1:N
        #=
        w1[i] = imag(expa/(-im*ω*σ)^(i))
        w2[i] = imag(expb/(-im*ω*σ)^(i))
        =#
        w1[i] = imag(expa*rp)
        w2[i] = imag(expb*rp)
        rp *= r
    end
    asymptotic_F!(a/σ, containers)
    Ia = dot(w1, aux) 
    asymptotic_F!(b/σ, containers)
    Ib = dot(w2, aux) 
    return Ia - Ib
end

function find_integration_interval(ϵ, a, b, σ, ω, containers::AsymptoticContainers)
    if asymptotic_error(a, b, σ, ω, containers) < ϵ return a end
    if a == σ a += (b-a)*0.01 end
    f = let b=b, σ=σ, ω=ω, containers=containers, ϵ=ϵ
        x -> asymptotic_error(x, b, σ, ω, containers) - ϵ
    end
    problem = ZeroProblem(f, (a, b))
    return solve(problem, Roots.Brent(), xatol=1e-0)
end
