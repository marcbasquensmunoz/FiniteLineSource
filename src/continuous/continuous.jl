@with_kw struct PrecomputationCoeffs{T <: Number}
    A::Vector{T}
    B::Vector{T}
    D::Vector{T}

    M::Array{T, 3}
    N::Matrix{T}

    fk
end

function legendre_coeffs(ff, x, w)
    # x should be in the range [-1, 1]
    # ff are the values of f in the range [a, b]
    n = length(x)
    c = zeros(n)
    for k in 0:n-1
        c[k+1] = (2k+1)/2 * sum([w[l]*Pl(x[l], k)*ff[l] for l in 1:n])
    end
    return c
end

function advance_1!(f, qk, A, B)
    N = length(qk)
    for k in 1:N
        f[:, k] = A .* f[:, k] + B .* qk[k]
    end
end
function advance_2!(f, qk, C)
    K = length(qk)
    for k in 1:K
        f[:, k] += C .* qk[k]
    end
end
function integral!(Tk; fk, qk, M, N, C)
    K = size(M)[2]
    n = size(M)[3]
    mult = reshape(reshape(M, K^2, n)*fk, K, K, K)
    for i in 1:K
        Tk[i] = C .* (sum([mult[i, j, j] for j in 1:K]) + (N*qk)[i]) 
    end
end


fR(r; a, m, c) = imag(im^a * exp(im*r*c)) / r^(3/2) * besselj(a+1/2, r*m)
function r̃(z1, z2, geometry::SegmentToSegment, rb) 
    @unpack D1, D2, H1, H2, σ = geometry
    sqrt(σ^2 + (H1/2*(z1+1) + D1 - H2/2*(z2+1) - D2)^2) / rb
end
function int_f_matrix_1(geometry::SegmentToSegment; n, K, m, c, rb, x, w)
    @unpack D1, D2, H1, H2, σ = geometry

    nk = length(x)
    # M_lka
    M = zeros(K+1, K+1, n+1)

    RR = zeros(nk, nk)

    for a in 0:n

        for i in 1:nk
            for j in 1:nk
                RR[i, j] = fR(r̃(x[i], x[j], geometry, rb), a=a, m=m, c=c)
            end
        end

        for k in 0:K
            for l in 0:K
                M[l+1, k+1, a+1] = (2l+1)/2 * sum([w[i]*w[j] * Pl(x[i], k) * Pl(x[j], l) * RR[i, j] for i in 1:nk, j in 1:nk])
            end
        end
    end
    M .* (H1/2)
end

function int_f_matrix_2(;x, w, m)
    n = length(x)
    sqrt(m*π/2) .* [(2a+1) * w[b] * Pl(x[b], a) for a in 0:n-1, b in 1:n]
end

function int_q_matrix(geometry::SegmentToSegment; K, rb, xx, ww)
    @unpack D1, D2, H1, H2, σ = geometry
    xN = length(xx)
    π/2 * H1/2 .* [(2l+1)/2 * sum([ww[i]*ww[j] * Pl(xx[i], k) * Pl(xx[j], l) / r̃(xx[i], xx[j], geometry, rb) for i in 1:xN, j in 1:xN]) for l in 0:K, k in 0:K]
end

function precompute_matrices(geometry::SegmentToSegment, params::Constants; n=200, K=10, nk=100)
    @unpack Δt, α, rb, b = params
    Δt̃ = α*Δt/rb^2
    a = 0.
    m = (b-a)/2
    c = (b+a)/2

    x, w = gausslegendre(n+1)
    ζ = @. m*x + c

    A = @. exp(-Δt̃ * ζ^2)
    B = @. -exp(-Δt̃ * ζ^2) / ζ
    D = @. 1/ζ

    xx, ww = gausslegendre(nk+1)

    R = int_f_matrix_1(geometry, n=n, K=K, m=m, c=c, rb=rb, x=xx, w=ww)
    P = int_f_matrix_2(x=x, w=w, m=m)
    M = reshape(reshape(R, (K+1)^2, n+1)*P, K+1, K+1, n+1)

    N = int_q_matrix(geometry, K=K, rb=rb, xx=xx, ww=ww)

    fk = zeros(n+1, K+1)

    PrecomputationCoeffs(A=A, B=B, D=D, M=M, N=N, fk=fk)
end

function compute_coefficients_through_history(Tk, qk, #=geometry::SegmentToSegment,=# params::Constants, precomp::PrecomputationCoeffs)
    @unpack A, B, D, M, N, fk = precomp
    @unpack kg, rb = params
    C = 1 / (2 * π^2 * kg * rb)
    T = size(Tk, 2)
    for i in 1:T
        advance_1!(fk, qk, A, B)
        integral!(@view(Tk[:, i]), fk=fk, qk=qk, M=M, N=N, C=C)
        advance_2!(fk, qk, D)
    end
end