
function create_error_plot(ϵ_range, r̃_range, res_step, res_synth; xlabel=:sigma, title)
    CairoMakie.activate!()
    fig = Figure()
    xlabel_text = xlabel == :sigma ? L"\log_{10} \ \tilde{\sigma}" : L"\log_{10} \ \tilde{r}"
    ax = Axis(fig[1, 1], xlabel = xlabel_text, ylabel = L"\log_{10} \ \Vert \varepsilon  \Vert_{\infty}", title = title)

    color_range = Makie.wong_colors()

    for (i, ϵ) in enumerate(ϵ_range)
        lines!(ax, log10.(r̃_range), log10.(res_step[:, i]), color=color_range[i])
        lines!(ax, log10.(r̃_range), log10.(res_synth[:, i]), color=color_range[i], linestyle = :dash)
    end
    Makie.xlims!(ax, minimum(log10.(r̃_range)), maximum(log10.(r̃_range)))
    Makie.ylims!(ax, -16., 0.)

    legend_colors = [PolyElement(color = color, strokecolor = :transparent) for color in color_range]
    legend_markers = [LineElement(color = :black), LineElement(color = :black, linestyle = :dash)]
    error_strings = [L"10^{%$(Int(log10(err)))}" for err in ϵ_range]
    legend = Legend(fig, [legend_colors, legend_markers], [error_strings, ["Unit step", "Synthetic"]], [L"\epsilon", "Load applied"], tellheight = false, tellwidth = true)
    legend.titleposition = :top
    legend.orientation = :vertical

    fig[1,2] = legend
    fig
end