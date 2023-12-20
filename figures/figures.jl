using FiniteLineSource
using Plots

r = 0.5
Δt = 3600.
n = 5000
a, b = 0., 5.

α, rb, kg = 10^-6, 0.1, 3.
Δt̃ = α*Δt/rb^2

# More complex load (intends to be realistic)
H = 8760
D = 24
t = 1:5H
q = [20*sin(2*π/H*i) + 5*sin(2*π/D*i) + 5 for i in t]

# Load
plot(t, q, label="")
xlabel!("t (h)")
ylabel!("q' (W/m)")

function plot_F(t)
    x = a:(b-a)/(n-1):b
    fx = zeros(n)
    for qt in q[1:t]
        FiniteLineSource.evolve!(fx, x, Δt̃, qt)
    end
    plot(x, fx, title="$t h", label="")
end
plot1 = plot_F(1)
plot2 = plot_F(10)
plot3 = plot_F(100)
plot4 = plot_F(1000)
plot5 = plot_F(10000)
plot6 = plot_F(40000)
#Plot of F at various times
plot(plot1, plot2, plot3, plot4, plot5, plot6, layout=(3,2))


# Plot T as a function of the distance
r = 0:0.1:50
Tr = [FiniteLineSource.compute_series(q, r=i, Δt=3600., n=400)[length(q)] for i in r]
plot(r ./ rb, Tr, label="",  ylims=(-0.03,0.5))
xlabel!("r̃")
ylabel!("T - T₀ (K)")

r = 0:1:50
function compute_rel_err_r(q, nv, dv)
    Tr = [FiniteLineSource.compute(q, r=i, Δt=3600., nv=nv, dv=dv) for i in r]
    Cr = [FiniteLineSource.convolve_step(q, Δt=Δt, r=i)[length(q)] for i in r]
    return @. log(abs((Tr - Cr) / Cr))
end
rel_err_r_1 = compute_rel_err_r(q, Int64.(100*ones(6)), [0., 0.1, 0.2, 0.5, 1., 5., 10.])
rel_err_r_2 = compute_rel_err_r(q, Int64.(200*ones(6)), [0., 0.1, 0.2, 0.5, 1., 5., 10.])
rel_err_r_3 = compute_rel_err_r(q, Int64.(400*ones(6)), [0., 0.1, 0.2, 0.5, 1., 5., 10.])
plot(r ./ rb, rel_err_r_1, label="n=100", ylims=(-35,0))
plot!(r ./ rb, rel_err_r_2, label="n=200")
plot!(r ./ rb, rel_err_r_3, label="n=400")
xlabel!("r̃")
ylabel!("log(ϵ)")


# Plot of temporal series of T
function compute_rel_err(q, r, b=10.)
    Tt = FiniteLineSource.compute_series(q, Δt=Δt, r=r, n=400, α=α, kg=kg, b=b)
    Ct = FiniteLineSource.convolve_step(q, Δt=Δt, r=r)
    return @. log(abs((Tt - Ct) / Ct))
end
rel_err_1 = compute_rel_err(q, 1)
rel_err_2 = compute_rel_err(q, 5)
rel_err_3 = compute_rel_err(q, 10, 7)
tind = 1:50:length(t)
plot(t[tind], rel_err_1[tind], title="n=400", label="r̃ = 10")
plot!(t[tind], rel_err_2[tind], label="r̃ = 50")
plot!(t[tind], rel_err_3[tind], label="r̃ = 100")
xlabel!("t (h)")
ylabel!("log(ϵ)")


# Plot of accuracy depending on n
N = 1:5:500
function compute_rel_err_n(q, r)
    res = FiniteLineSource.convolve_step(q, Δt=Δt, r=r)[length(q)]
    Tn = [FiniteLineSource.compute(q, Δt=Δt, r=r, α=α, kg=kg, nv=Int64.(n*ones(6)), dv=[0., 0.1, 0.2, 0.5, 1., 5., 10.]) for n in N]
    return @. log(abs((Tn - res) / res))
end
rel_err_n_1 = compute_rel_err_n(q, 1)
rel_err_n_2 = compute_rel_err_n(q, 5)
rel_err_n_3 = compute_rel_err_n(q, 10)
plot(N, rel_err_n_1, label="r̃ = 10")
plot!(N, rel_err_n_2, label="r̃ = 50")
plot!(N, rel_err_n_3, label="r̃ = 100")
xlabel!("n")
ylabel!("log(ϵ)")