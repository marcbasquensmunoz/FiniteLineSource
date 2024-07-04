using FiniteLineSource
using WGLMakie
using CairoMakie


function produce_plot_fls(σ, line_points, plot)        
    q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]
    I = zeros(length(q))

    setup = SegmentToPoint(D=0., H=10., σ=σ, z=5.)
    params = Constants(Δt = 3600., line_points=line_points, line_limits=[0., 0.5, 1.])#[0., 0.3, 0.5, 0.7, 1.])

    precomp = precompute_parameters(setup, params=params)
    compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    C = convolve_step(q, setup, params=params)
    minerr = 5*10^-17

    abs_error = abs.(C-I) 
    abs_error[abs_error .< minerr] .= minerr
    error = log10.(abs_error)

    density!(plot, error, direction = :y)
end

####

rs = [0.1, 1, 10, 100]
points = [1, 2, 3, 4, 6]
#discretization = [20, 30, 30, 20]
discretization = [10, 10]
ff = Figure(size = (1000,1000))
grid = ff[1, 1] = GridLayout()
axis = Matrix{Axis}(undef, length(rs), length(points))

xlimits = (0, 4.2)
limits = (-17, 0)

for (i, r) in enumerate(rs)
    for (j, p) in enumerate(points)
        axis[i, j] = Axis(grid[i, j], limits=(xlimits, limits), ylabel=j==1 ? L"\mathbf{\tilde{\sigma}=%$(Int(r/0.1))}" : "", ylabelrotation = 0)
        real_points = produce_plot_fls(r, p .* discretization, axis[i, j])    
        axis[i, j].title = i==1 ? L"\mathbf{\sum N = %$(Int(p * sum(discretization) + length(discretization)))}" : ""
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
save("figures/error_stp.pdf", ff)