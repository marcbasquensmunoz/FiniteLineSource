using FiniteLineSource
using FiniteLineSource: convolve_step
using WGLMakie
using CairoMakie
using LaTeXStrings

function produce_plot(r, n, plot)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = PointToPoint(r=r)

    params = Constants(Δt = 3600.)

    n_tot = [0]
    precomp = precompute_parameters(setup, params=params, n=n, n_tot=n_tot)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)
    
    C = convolve_step(q, Δt = params.Δt, r = setup.r)
    minerr = 5*10^-17

    abs_error = abs.(C-I) 
    abs_error[abs_error .< minerr] .= minerr
    error = log10.(abs_error)

    density!(plot, error, direction = :y)
    return n_tot[1]
end

####

rs = [0.1, 1, 10, 100]
points = [4, 8, 12, 16, 20]
ff = Figure(size = (1000,1000))
grid = ff[1, 1] = GridLayout()
axis = Matrix{Axis}(undef, length(rs), length(points))

xlimits = (0, 4.2)
limits = (-17, 0)

for (i, r) in enumerate(rs)
    for (j, p) in enumerate(points)
        axis[i, j] = Axis(grid[i, j], limits=(xlimits, limits), ylabel=j==1 ? L"\mathbf{\tilde{r}=%$(Int(r/0.1))}" : "", ylabelrotation = 0)
        real_points = produce_plot(r, p, axis[i, j])    
        axis[i, j].title = i==1 ? L"\mathbf{\sum n = %$(real_points)}" : ""
    end     
end
for i in 1:length(rs)-1
    for j in 1:length(points)-1
        linkxaxes!(axis[i, j], axis[i+1, j])
        linkyaxes!(axis[i, j], axis[i, j+1])
    end
end

for i in 1:length(rs)-1
    for j in 1:length(points)
        hidexdecorations!(axis[i, j], grid = false)
    end
end

for j in 2:length(points)
    for i in 1:length(rs)
        hideydecorations!(axis[i, j], grid = false)
    end
end

colgap!(grid, 5)
rowgap!(grid, 5)

ff

CairoMakie.activate!()
save("figures/error_ps.pdf", ff)
