using Plots
using FiniteLineSource
using FiniteLineSource: point_step_response, step_response
using QuadGK
using FastGaussQuadrature
using LinearAlgebra

Ds = 0.
Hs = 10.

Dt = 0.
Ht = 10.

σ = 0.1
t = 3600.#5*24*3600.

mf = 1.
cf = 4000.

α = 1e-6
kg = 3.

Tfinput = 300.

Q(z) = 0.05 * (Hs^2/4 - (z-Hs/2)^2) + 1
plot(Ds:0.01:Ds+Hs, Q.(Ds:0.01:Ds+Hs))

n = 6   # Segments
m = 100  # Points per segment

params = Constants()
A = [-0.00396062862871448 -0.0002694126315114766; 0.0002694126315114766 0.00396062862871448]

# Shape of weight function
L = 10
xx = 0:0.01:L
f(z) =  sum([1 -1]*exp(A*(L-z))*A*[1;1]) / sum([1 -1]*exp(A*L)*[1;1])
plot(xx, f.(xx))
##

I1 = quadgk(z1 -> quadgk(z2 -> f(z2)*point_step_response(t, sqrt(σ^2 + (z1-z2)^2), α, kg), 0, L)[1], 0, L)
I2 = quadgk(z1 -> quadgk(z2 ->  1/L *point_step_response(t, sqrt(σ^2 + (z1-z2)^2), α, kg), 0, L)[1], 0, L)
##

zz = zeros(n*m)
ww = zeros(n*m)
x, w = gausslegendre(m)

for i in 1:n
    t1 = Dt + (i-1)*Ht/n
    t2 = Dt + i*Ht/n
    zz[m*(i-1)+1:m*i] = (t2-t1)/2 .* x .+ (t2+t1)/2 
    ww[m*(i-1)+1:m*i] = w .* (t2-t1)/2
end

Ez = map(z -> exp(A*z), zz)

Tb_real(z) = quadgk(ζ -> Q(ζ) * point_step_response(t, sqrt(σ^2 + (z-ζ)^2), α, kg), Ds, Ds+Hs, atol=1e-16)[1] 

Tbint = [quadgk(z -> exp(A*(Dt + i*Ht/n-z))*A*[1,1] .* Tb_real(z), Dt + (i-1)*Ht/n, Dt + i*Ht/n, atol=1e-16)[1] for i in 1:n]

function Tb_new_s(i)
    t1 = Dt + (i-1)*Ht/n
    t2 = Dt + i*Ht/n
    tt(i) = i == 0 ? 0. : Dt + i*Ht/n

    e = sum([1 -1] * exp(A*(t2-t1)) * [1, 1])
    tot = 0
    for j in 1:(i-1)
        tot += sum([1 -1] * (exp(A*(tt(i) - tt(j))) - exp(A*(tt(i-1) - tt(j)))) * (Tb_new_val[j] * (I - exp(A*(tt(j) - tt(j-1)))) * [1, 1] +  Tbint[j]))
    end
    tot += sum([1 -1]*Tbint[i])
    
    tot/e
end
Tb_new_val = zeros(n)
for i in 1:n
    Tb_new_val[i] = Tb_new_s(i)
end
function Tb_new(z)
    for i in 1:n
        z1 = Dt + (i-1)*Ht/n
        z2 = Dt + i*Ht/n
        if z1 <= z && z <= z2 
            return Tb_new_val[i]
        end
    end
    0.
end


function Tb_new_simple_s(i)
    t1 = Dt + (i-1)*Ht/n
    t2 = Dt + i*Ht/n

    e = sum([1 -1] * exp(A*(t2-t1)) * [1, 1])
    sum([1 -1] * quadgk(z -> exp(A*(t2-z)) .* Tb_real(z), t1, t2, atol=1e-16)[1] * A * [1, 1]) / e
end
Tb_new_simple_val = [Tb_new_simple_s(i) for i in 1:n]
function Tb_new_simple(z)
    for i in 1:n
        z1 = Dt + (i-1)*Ht/n
        z2 = Dt + i*Ht/n
        if z1 <= z && z <= z2 
            return Tb_new_simple_val[i]
        end
    end
    0.
end

function Tb_sts(i, j)
    s1 = Ds + (i-1)*Hs/n
    s2 = Ds + i*Hs/n

    t1 = Dt + (j-1)*Ht/n
    t2 = Dt + j*Ht/n

    q = quadgk(z -> Q(z), s1, s2)[1] / (s2 - s1)
    q * step_response(t, SegmentToSegment(D1=s1, H1=s2-s1, D2=t1, H2=t2-t1, σ=σ), params)
end

Tb_sts_val = [Tb_sts(i, j) for i in 1:n, j in 1:n]

function Tb_mean(z)
    for j in 1:n
        z1 = Dt + (j-1)*Ht/n
        z2 = Dt + j*Ht/n
        if z1 <= z && z <= z2 
            return sum(Tb_sts_val[:, j])
        end
    end
    0.
end


Tbz_real = @. Tb_real(zz)
Tbz_sts = @. Tb_mean(zz)
Tbz_new = @. Tb_new(zz)
Tbz_simple_new = @. Tb_new_simple(zz)

plot(zz, Tbz_real)
plot!(zz, Tbz_sts)
plot!(zz, Tbz_new) 
plot!(zz, Tbz_simple_new) 

plot(zz, Tbz_new - Tbz_sts)

EH = exp(A*Ht) 
Ein = EH[2, 1] - EH[1, 1]
Eout = EH[1, 2] - EH[2, 2]


function compute_Q(Tbz)
    Tfoutput = (Ein*Tfinput + sum([1 -1] * sum(ww .* Ez[end:-1:1] .* Tbz[:]) * A * [1, 1]))/Eout
    Tf0 = [Tfinput, Tfoutput]
    Tfs = map(i -> Ez[i*m]*Tf0 - sum([sum(ww[(j-1)*m+1:j*m] .* Ez[m*(i-j+1):-1:m*(i-j)+1].* Tbz[(j-1)*m+1:j*m])*A*[1,1] for j in 1:i]), 1:n)
    pushfirst!(Tfs, Tf0)

    Qs = [sum(mf*cf*[1 -1]*(Tfs[i] - Tfs[i-1])) for i in 2:n+1]
    return Qs
end

Qs_real = compute_Q(Tbz_real)
Qs_sts = compute_Q(Tbz_sts)
Qs_new = compute_Q(Tbz_new)
Qs_simple_new = compute_Q(Tbz_new)


Qs_real-Qs_sts
Qs_real-Qs_new
Qs_sts-Qs_new

Qs_new - Qs_simple_new

@show sum(abs.(Qs_real-Qs_new))
@show sum((Qs_real-Qs_new))
@show sum(abs.(Qs_real-Qs_sts))
@show sum((Qs_real-Qs_sts))


function error(i)
    tt(j) = j == 0 ? 0. : Dt + j*Ht/n
    t1 = tt(i-1)
    t2 = tt(i)

    e = sum([1 -1] * exp(A*(t2-t1)) * [1, 1])

    I1 = quadgk(z -> exp(-A*z) * A .* Tb_real(z), 0, t1, atol=1e-16)[1]
    I2 = i > 1 ? sum([Tb_new_val[j] .* (exp(-A*tt(j)) - exp(-A*tt(j-1))) for j in 1:i-1]) : zeros(2, 2)
    sum([1 -1] * (exp(A*t2) - exp(A*t1)) * ( I1 + I2 ) * [1, 1]) / e
end

error.(1:n)

#### h of new 
i=2
tt(i) = i == 0 ? 0. : Dt + i*Ht/n
params_mean = MeanSegToSegEvParams(D1=Ds, H1=Ds+Hs, D2=tt(i-1), H2=tt(i)-tt(i-1), σ=σ)
params = InternalSegToSegEvParams(D1=Ds, H1=Ds+Hs, D2=tt(i-1), H2=tt(i)-tt(i-1), σ=σ, A=A)

h_mean, r_min, r_max = mean_sts_evaluation(params_mean)
h_new, _, _ = mean_internal_evaluation(params)

rr = r_min:0.01:r_max

plot(rr, h_mean.(rr))
plot!(rr, h_new.(rr))

diff_h = h_mean.(rr) -  h_new.(rr)
sum(abs.(diff_h[2:end]))

plot(rr, diff_h)