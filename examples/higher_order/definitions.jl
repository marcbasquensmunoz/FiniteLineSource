using Parameters 
using FastGaussQuadrature
using QuadGK
using LinearAlgebra
using FiniteLineSource

indices(i::Int, N::Int) = (div(i-1, N)+1, (i-1)%N+1)
struct LagrangeBasis{T <: Number}
    coefs::Matrix{T}
    x::Vector{T}
    w::Vector{T}
    N::Int
end
function LagrangeBasis(x, w)
    N = length(x)
    Ξ = [x[k]^m for k in 1:N, m in 0:N-1]
    δ = Diagonal(ones(N))
    # Is this faster than solving N systems of equations?
    # A = inv(Ξ)*δ
    A = zeros(N, N)
    for i in 1:N
        @views A[:, i] = Ξ \ δ[:, i]
    end
    LagrangeBasis(A, x, w, N)
end
function LagrangeBasis(N::Int) 
    x, w = gausslegendre(N)
    LagrangeBasis(x, w)
end
function ψ(ξ, k::Int, b::LagrangeBasis) 
    @unpack coefs, N = b
    res = 0.
    for i in N:-1:1
        res *= ξ
        res += coefs[i, k]
    end
    return res
end

@with_kw struct InternalModelParams{T <: Number} @deftype T
    mf
    cpf
    R11
    R22
    R12
    Rp = (R11+R22)/(R11*R22)

    β1 = 1/(mf*cpf*R11)
    β2 = 1/(mf*cpf*R22)
    β12 = 1/(mf*cpf*R12)
    
    β = (β2-β1)/2
    γ = sqrt((β2+β1)^2 / 4 + (β2+β1)*β12)
    δ = 1/γ * (β12 + (β2+β1)/2)
end

struct BoreholeDiscretization{T <: Number}
    D::T
    H::T
    segments::Vector{T}
    L::Vector{T}
    S::Int
end
BoreholeDiscretization(D, H, segments) = BoreholeDiscretization(D, H, segments, (segments[2:end] - segments[1:end-1]) .* H/2, length(segments) - 1)

s(ξ, p) = 50. * (ξ+1) #H/2 * (ξ+1) + D
f1(ξ, p::InternalModelParams) = exp(p.β*s(ξ, p)) * (cosh(p.γ*s(ξ, p)) - p.δ*sinh(p.γ*s(ξ, p)))
f2(ξ, p::InternalModelParams) = exp(p.β*s(ξ, p)) * p.β12/p.γ * sinh(p.γ*s(ξ, p))
f3(ξ, p::InternalModelParams) = exp(p.β*s(ξ, p)) * (cosh(p.γ*s(ξ, p)) + p.δ*sinh(p.γ*s(ξ, p)))
f4(ξ, p::InternalModelParams) = exp(p.β*s(ξ, p)) * (p.β1*cosh(p.γ*s(ξ, p)) - (p.δ*p.β1 + p.β2*p.β12/p.γ)*sinh(p.γ*s(ξ, p)))
f5(ξ, p::InternalModelParams) = exp(p.β*s(ξ, p)) * (p.β2*cosh(p.γ*s(ξ, p)) + (p.δ*p.β2 + p.β1*p.β12/p.γ)*sinh(p.γ*s(ξ, p)))

ξseg(ξ′, u, segments) = (segments[u+1] - segments[u]) / 2 * ξ′ +  (segments[u+1] + segments[u]) / 2

function evaluate(ξ, coefs, bh_disc::BoreholeDiscretization, basis::LagrangeBasis) 
    @unpack N = basis
    @unpack segments, S = bh_disc
    for k in 1:S
        if segments[k] <= ξ <= segments[k+1]
            return sum([coefs[(k-1)*N+i] * ψ(2*(ξ-segments[k])/(segments[k+1]-segments[k]) - 1, i, basis) for i in 1:N])
        end
    end
end