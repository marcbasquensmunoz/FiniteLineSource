using FiniteLineSource
using WGLMakie
using CairoMakie
using Cubature


function produce_plot_sts(line_points, plot, q, setup, C)        
    I = zeros(length(q))

    params = Constants(Δt = 3600., line_points=line_points, line_limits=[0., 0.3, 0.5, 0.7, 1.], )
    precomp = precompute_parameters(setup, params=params)
    @time compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

    minerr = 5*10^-17

    abs_error = @. abs.(C-10I) 
    @show maximum(abs.(C-I))
    @show maximum(abs.(C- 10 .* I))
    abs_error[abs_error .< minerr] .= minerr
    error = log10.(abs_error)

    density!(plot, error, direction = :y)
end

####

rs = [0.1, 1., 10., 100.]
points = [1, 2, 3, 4, 20]
discretization = [4, 6, 6, 4]
ff = Figure(size = (1000,1000))
grid = ff[1, 1] = GridLayout()
axis = Matrix{Axis}(undef, length(rs), length(points))

xlimits = (0, 5.2)
limits = (-17, 0)
C = zeros(8760*20, length(rs))

q = [20*sin(2π*i/8760) + 5*sin(2π*i/24) + 5. for i=1:8760*20]

for (i, r) in enumerate(rs)
    setup = SegmentToSegmentOld(D1=0., H1=10., D2=0., H2=10., σ=r)

    params = Constants(Δt = 3600.)
    C = @time convolve_step(q, setup, params=params)

    for (j, p) in enumerate(points)
        @show r, p
        axis[i, j] = Axis(grid[i, j], limits=(xlimits, limits), ylabel=j==1 ? L"\mathbf{\tilde{\sigma}=%$(Int(r/0.1))}" : "", ylabelrotation = 0)
        real_points = produce_plot_sts(p .* discretization, axis[i, j], q, setup, C)    
        axis[i, j].title = i==1 ? L"\mathbf{\sum N = %$(Int(p * sum(discretization) + length(discretization)))}" : ""
        if i == 4 && j == 3
            axis[i, j].xlabel = "Density of the error "
        end
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
save("figures/error_sts.pdf", ff)