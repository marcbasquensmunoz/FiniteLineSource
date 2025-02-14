using Parameters 
using FastGaussQuadrature
using QuadGK
using LinearAlgebra
using FiniteLineSource
using StaticArrays

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
    R11Δ
    R22Δ
    R12Δ
    Rp = (R11Δ+R22Δ)/(R11Δ*R22Δ)
    Rb = (R11Δ*R22Δ)/(R11Δ+R22Δ)

    β1 = 1/(mf*cpf*R11Δ)
    β2 = 1/(mf*cpf*R22Δ)
    β12 = 1/(mf*cpf*R12Δ)
    
    β = (β2-β1)/2
    γ = sqrt((β2+β1)^2 / 4 + (β2+β1)*β12)
    δ = 1/γ * (β12 + (β2+β1)/2)
end

struct BoreholeDiscretization{T <: Number}
    #η::Vector{T}
    #p::Vector{Vector{T}}
    A::Matrix{T}
    Nd::Int
    segments::Vector{T}
    S::Int
end
function BoreholeDiscretization(ξ, p, segments = ξ) 
    Nd = 2
    P = [map(x->x[i], p) for i in 1:3]
    Ξ = [ξ[k]^m for k in 1:Nd, m in 0:Nd-1]

    A = zeros(Nd, 3)
    for i in 1:3
        A[:, i] = Ξ \ P[i]
    end
    BoreholeDiscretization(A, Nd, segments, length(segments) - 1)
end
function path(ξ, A, Nd)
    x = 0.
    y = 0.
    z = 0.
    for i in Nd-1:-1:0
        x *= ξ
        y *= ξ
        z *= ξ
        x += A[i+1, 1]
        y += A[i+1, 2]
        z += A[i+1, 3]
    end
    return @SVector [x, y, z]
end
path(ξ, d::BoreholeDiscretization) = path(ξ, d.A, d.Nd)
function path(ξ, segment::Int, d::BoreholeDiscretization) 
    @unpack segments = d
    m = (segments[segment+1] - segments[segment])/2
    c = (segments[segment+1] + segments[segment])/2
    path(m*ξ+c , d)
end
function dpath(ξ, d::BoreholeDiscretization)
    @unpack A, Nd = d
    x = 0.
    y = 0.
    z = 0.
    for i in Nd-1:-1:1
        x *= ξ
        y *= ξ
        z *= ξ
        x += convert(eltype(A), i) * A[i+1, 1]
        y += convert(eltype(A), i) * A[i+1, 2]
        z += convert(eltype(A), i) * A[i+1, 3]
    end
    return @SVector [x, y, z]
end

s(ξ, d::BoreholeDiscretization) = quadgk(η -> norm(dpath(η, d)), -1., ξ)[1]
function normJ(η, segment::Int, d::BoreholeDiscretization)
    @unpack segments = d
    m = (segments[segment+1] - segments[segment])/2
    c = (segments[segment+1] + segments[segment])/2
    m * norm(dpath(m*η + c, d))
end

function f1(ξ, p::InternalModelParams, d::BoreholeDiscretization) 
    @unpack β, δ, γ = p
    sl = s(ξ, d)
    exp(β*sl) * (cosh(γ*sl) - δ*sinh(γ*sl))
end
function f2(ξ, p::InternalModelParams, d::BoreholeDiscretization) 
    @unpack β, β12, γ = p
    sl = s(ξ, d)
    exp(β*sl) * β12/γ * sinh(γ*sl)
end
function f3(ξ, p::InternalModelParams, d::BoreholeDiscretization) 
    @unpack β, γ, δ = p
    sl = s(ξ, d)
    exp(β*sl) * (cosh(γ*sl) + δ*sinh(γ*sl))
end
function f4(ξ1, ξ2, p::InternalModelParams, d::BoreholeDiscretization) 
    @unpack β, β1, β2, β12, γ, δ = p
    sl = s(ξ1, d) - s(ξ2, d)
    exp(β*sl) * (β1*cosh(γ*sl) - (δ*β1 + β2*β12/γ)*sinh(γ*sl))
end
function f5(ξ1, ξ2, p::InternalModelParams, d::BoreholeDiscretization) 
    @unpack β, β1, β2, β12, γ, δ = p
    sl = s(ξ1, d) - s(ξ2, d)
    exp(β*sl) * (β2*cosh(γ*sl) + (δ*β2 + β1*β12/γ)*sinh(γ*sl))
end

function f1p2(ξ, p::InternalModelParams, d::BoreholeDiscretization; C1 = 1., C2 = 1.) 
    @unpack β, β12, γ, δ = p
    sl = s(ξ, d)
    exp(β*sl) * ( C1 * cosh(γ*sl) + (C2 * β12/γ - C1 * δ ) * sinh(γ*sl))
end
function f2p3(ξ, p::InternalModelParams, d::BoreholeDiscretization; C2 = 1., C3 = 1.) 
    @unpack β, β12, γ, δ = p
    sl = s(ξ, d)
    exp(β*sl) * ( C3 * cosh(γ*sl) + (C2 * β12/γ + C3 * δ ) * sinh(γ*sl))
end
function f4p5(ξ1, ξ2, p::InternalModelParams, d::BoreholeDiscretization; C4 = 1., C5 = 1.) 
    @unpack β, β1, β2, β12, γ, δ = p
    sl = s(ξ1, d) - s(ξ2, d)
    exp(β*sl) * ((C4*β1+C5*β2)*cosh(γ*sl) + (δ*(C5*β2-C4*β1) + β12/γ*(C5*β1-C4*β2) )*sinh(γ*sl))
end

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

function ξ_discretization(bh::BoreholeDiscretization, basis::LagrangeBasis)
    @unpack segments, S = bh
    @unpack x = basis
    reduce(vcat, [ξseg.(x, Ref(u), Ref(segments)) for u in 1:S])
end

function integrate_bh(f, bh::BoreholeDiscretization, basis::LagrangeBasis)
    gQ = zeros(length(f))
    g_Q!(gQ, bh, basis)
    dot(gQ, f)    
end
