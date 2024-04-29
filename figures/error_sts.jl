using FiniteLineSource
using WGLMakie
using CairoMakie


function produce_plot_fls(r, line_points, full, dist, year)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = SegmentToSegmentOld(D1=0., H1=20., D2=0., H2=20., r=r)
    params = Constants(Δt = 3600., line_points=line_points, line_limits=[0., 0.3, 0.5, 0.7, 1.])

    precomp = precompute_parameters(setup, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, setup, params=params)
    plot_start_year = 5*8760
    plot_step = 24
    plot_end_year = 6*8760
    minerr = 5*10^-17

    abs_error = abs.(C-I) 
    t = (1:plot_step:length(q)) ./ 8760
    t_year = (0:8760) ./ 8760 * 12

    abs_error[abs_error .< minerr] .= minerr
    error = log10.(abs_error)

    lines!(full, t, error[1:plot_step:end], label="∑N=$(Int(sum(line_points)))")
    density!(dist, error, direction = :y)
    lines!(year, t_year, error[plot_start_year:plot_end_year])
end

rs = [0.1, 1, 10, 100]
ff = Figure(size = (1000,1000))
grid = ff[1, 1] = GridLayout()
axis = Matrix{Axis}(undef, length(rs), 3)

limits = (-17, 0)

for (i, r) in enumerate(rs)
    axis[i, 1] = Axis(grid[i, 1], limits=(nothing, limits), xlabel=i==length(rs) ? "t (years)" : "",  title=i==1 ? L"\mathbf{\log_{10}(\epsilon(t))}" : "", ylabel=L"\mathbf{\tilde{r}=%$(Int(r/0.1))}", ylabelrotation = 0)
    axis[i, 2] = Axis(grid[i, 2], limits=(nothing, limits), title=i==1 ? L"\mathbf{\log_{10}(\epsilon)} \text{ distribution}" : "")
    axis[i, 3] = Axis(grid[i, 3], limits=(nothing, limits), xlabel=i==length(rs) ? "t (months)" : "",  title=i==1 ? L"\mathbf{\log_{10}(\epsilon(t))} \text{ through a cycle}" : "")

    linkyaxes!(axis[i, 1], axis[i, 2])
    linkyaxes!(axis[i, 2], axis[i, 3])
    hideydecorations!(axis[i, 2], grid = false)
    hideydecorations!(axis[i, 3], grid = false)

    produce_plot_fls(r, 0.5 .* [20, 30, 30, 20], axis[i, :]...)    
    produce_plot_fls(r, [20, 30, 30, 20], axis[i, :]...)       
    produce_plot_fls(r, 2 .* [20, 30, 30, 20], axis[i, :]...)       
    produce_plot_fls(r, 3 .* [20, 30, 30, 20], axis[i, :]...)     
    produce_plot_fls(r, 4 .* [20, 30, 30, 20], axis[i, :]...)          
end
for i in 1:length(rs)-1
    for j in 1:3
        linkxaxes!(axis[i, j], axis[i+1, j])
        hidexdecorations!(axis[i, j], grid = false)
    end
end
axis[length(rs), 3].xticks = 0:3:12

colsize!(grid, 2, Fixed(150))
colgap!(grid, 5)
rowgap!(grid, 5)
leg = Legend(grid[4, 4], axis[1, 1], linewidth=5.)

ff

CairoMakie.activate!()
save("figures/error_sts.pdf", ff)