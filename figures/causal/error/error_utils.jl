
function create_error_plot(ϵ_range, r_range, res_step, res_synth; xlabel=:sigma, title, rb)
    CairoMakie.activate!()
    fig = Figure()
    xlabel_text = xlabel == :sigma ? L"\tilde{\sigma}" : L"\tilde{r}"
    ax = Axis(fig[1, 1], xlabel = xlabel_text, ylabel = L"\log_{10} \ \Vert \varepsilon  \Vert_{\infty}", title = title)

    for (i, ϵ) in enumerate(ϵ_range)
        lines!(ax, r_range ./ rb, log10.(res_step[:, i]), label = L"\epsilon = 10^{%$(Int(log10(ϵ)))}")
        lines!(ax, r_range ./ rb, log10.(res_synth[:, i]), linestyle = :dash,  label = L"\epsilon = 10^{%$(Int(log10(ϵ)))}")
    end
    Makie.xlims!(ax, 0, maximum(r_range) / rb)
    Makie.ylims!(ax, -16., 0.)
    axislegend(""; position= :rt, backgroundcolor = (:grey90, 0.25));

    fig
end