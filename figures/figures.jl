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
Tr = [FiniteLineSource.compute(q, r=i, Δt=3600.)[length(q)] for i in r]
plot(r ./ rb, Tr, label="",  ylims=(-0.03,0.5))
xlabel!("r̃")
ylabel!("T - T₀ (K)")


# Plot of temporal series of T
Tt = FiniteLineSource.compute_series(q, Δt=Δt, r=1, n=300, α=α, kg=kg)
Ct = FiniteLineSource.convolve_step(q, Δt=Δt, r=1)
@. rel_err = log(abs((Tt - Ct) / Ct))
plot(t, rel_err, title="Log of relative error, n=300", label="")
xlabel!("t (h)")
ylabel!("log(ϵ)")


# Plot of accuracy depending on n
N = 1:500
res = FiniteLineSource.convolve_step(q, Δt=Δt, r=1)[length(q)]
Tn = [FiniteLineSource.compute_series(q, Δt=Δt, r=1, n=n, α=α, kg=kg)[length(q)] for n in N]
plot(N, log.(abs.(Tn .- res) ./ res), label="Relative error")
xlabel!("n")
ylabel!("log(ϵ)")