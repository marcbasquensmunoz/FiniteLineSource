using FiniteLineSource
using FiniteLineSource: PointSource
using Parameters
using LinearAlgebra
using SpecialFunctions
using Plots
using Statistics
using CairoMakie

function make_full_simulation(;Nt, r, ϵ, positions, constants, q)
    block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);

    Nb = size(block.K_min)[1]
    Ib = zeros(Nb, Nt)
    evolve!(Ib, q, block)
    return Ib
end

function make_interpolated_simulation_with_bound(;Nt, r, ϵ, positions, constants, q, windows = [], use_bound=false)
    block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);
    @unpack load_buffer, load_delays, qaux, expNin, expNout, ranges, expt, F, K_min, Kranges, HM, ζ, N = block
    @unpack Δt̃, Δt, α, kg, rb = constants

    Nb = size(K_min)[1]
    K = size(load_delays)[1]

    ϵ = ϵ / K

    for i in 1:K
        push!(windows, Int[])
    end

    I_int = zeros(Nb, Nt)
    N_start = 1
    Ns = 0.
    current_q = zeros(Nb)

    ref_source = 1
    ref_target = 2

    derf = [( erf(r/rb / sqrt(4*Nb*Δt̃)) - erf(r/rb / sqrt(4*(Nb+1)*Δt̃))) / (4*π*kg*r * Δt̃) for Nb in N]
    d2erf = [(erf(r/rb / sqrt(4*Nb*Δt̃)) + erf(r/rb / sqrt(4*(Nb+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(Nb+1)*Δt̃))) / (4*π*kg*r * Δt̃^2) for Nb in N]

    F_1 = zeros(size(F))

    q_in_buffer = zeros(K)
    q_out_buffer = zeros(K)

    q_in_der = zeros(2, K)
    q_out_der = zeros(2, K)

    dQ_in = zeros(K)
    dQ_out = zeros(K)

    Q_eff = zeros(K)

    bound_int = zeros(K)
    osc_int = zeros(K)
    C_disc_1 = zeros(K)
    C_disc_2 = zeros(K)
    B_tot = zeros(K)
    B_max = zeros(K)
    B = zeros(K)
    Ns = zeros(Int, K)

    N_start = N[1:K]

    I_k_past = zeros(K, Nb)
    I_k_present = zeros(K, Nb)

    for nt in 1:Nt
        for i in 1:Nb
            current_q[i] = load_buffer[i][1]
        end
        @views for (q, b) in zip(q[:,nt], load_buffer)
            push!(b, q)
        end

        @views @. F_1 = F
        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                q_in_buffer[i] = qin
                q_out_buffer[i] = qout
                push!(load_delays[i, j], qin)
                current_q[j] = qout
                @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
            end
            @views @. F[:, j] = expt * F[:, j] + qaux
        end
        
        if nt < N[1] || length(Kranges) == 0 continue end

        @views @. q_in_der[1, :] = q_in_der[2, :] 
        @views @. q_in_der[2, :] = q_in_buffer 

        @views @. q_out_der[1, :] = q_out_der[2, :] 
        @views @. q_out_der[2, :] = q_out_buffer 

        @views @. dQ_in = (q_in_der[2, :] - q_in_der[1, :]) / Δt̃
        @views @. dQ_out = (q_out_der[2, :] - q_out_der[1, :]) / Δt̃
        @views Q_in = q_in_der[2, :]
        @views Q_out = q_out_der[2, :]

        @. Q_eff = abs(Q_in)

        for i in 1:K
            osc_int[i] = dot(F[ranges[i], ref_source] .* (1 .- expt[ranges[i]]) .^ 2, HM[ranges[i], ref_target, ref_source]) / Δt̃^2
            C_disc_1[i] = derf[i] * dQ_in[i] - derf[i+1] * dQ_out[i]
            C_disc_2[i] = -d2erf[i] * Q_in[i] + d2erf[i+1] * Q_out[i]
            bound_int[i] = (log((N[i]+1)/(N[i+1]+1)) + 1/2 * log(N[i+1]*(N[i+1]+2)/(N[i]*(N[i]+2)))) / (2*π^2*kg*r*Δt̃^2) * Q_eff[i]
        end

        @. B_tot = (use_bound ? bound_int : osc_int) + C_disc_1 + C_disc_2
        @. B_max = max(abs(B_max), abs(B_tot))
        @. B = 1/Δt̃ * sqrt(8ϵ / B_max)
        @. Ns = [max(b == Inf ? 1 : Int(floor(b)), 1) for b in B]

        for i in 1:K
            if nt - N_start[i] >= Ns[i] || nt == Nt
                if nt > N[i] push!(windows[i], Ns[i]) end
                I_k_present[i, :] .= 0.
                for target in 1:Nb
                    for source in 1:Nb
                        if source == target continue end
                        @inbounds range = ranges[i]
                        @inbounds @views I_k_present[i, target] += dot(F[range, source], HM[range, target, source])
                    end
                    N_int = nt - N_start[i]
                    @views @. I_int[target, N_start[i]+1:nt] += I_k_past[i, target] + (I_k_present[i, target] - I_k_past[i, target]) * (1:N_int) / N_int
                end 
                @views I_k_past[i, :] .= I_k_present[i, :]

                N_start[i] = nt
                Ns[i] = 0
                B_max[i] = 0.
            end
        end
    end
    return I_int
end

ϵ = 1e-2
Δt = 3600.
Nt = 8760*200

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

r = 100.
positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
q = ones(2) * [0 + sin(i/24)  for i in 1:Nt]'
#q = 3 .* ones(2, Nt)

windows = Vector{Int64}[]

block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);
Ib = @time make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
I_int = @time make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows, use_bound=false)

maximum(abs.(Ib - I_int))
Plots.plot(log10.(abs.(Ib - I_int))[1,:])


Plots.plot(I_int[1,:])
Plots.plot!(Ib[1,:])

function get_max_error(q, ϵ, r; use_bound)
    Δt = 3600.
    Nt = size(q)[2]

    α = 1e-6
    kg = 3.
    rb = 0.1
    constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

    positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
    windows = Vector{Int64}[]

    Ib = make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
    I_int_bound = make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows, use_bound=use_bound)
    
    mean(abs.(Ib - I_int_bound)), windows
end

function generate_analysis_plot(q; use_bound=false, titles = [])
    rr = 10. .^ (0:0.2:2)
    ϵϵ = 10. .^ (-2:-2:-10)

    errors = zeros(length(rr), length(ϵϵ), length(q))
    windows = Array{Vector{Vector{Int}}}(undef, length(rr), length(ϵϵ), length(q))

    for (i, r) in enumerate(rr)
        for (j, ϵ) in enumerate(ϵϵ)
            for (k, load) in enumerate(q)
                @show k, r, ϵ
                error, w = get_max_error(load, ϵ, r, use_bound=use_bound)

                errors[i, j, k] = error
                windows[i, j, k] = w
            end
        end
    end
    errors[errors .< 1e-16] .= 1e-16

    max_K = 3

    axes = Matrix{Axis}(undef, max_K+1, length(q))
    fig = Figure(size=(length(q)*250, 900))

    for (i, load) in enumerate(q) 
        axes[1, i] = Axis(fig[1, i], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} \ \Vert \epsilon  \Vert_{\infty}")
        for j in 1:max_K
            axes[1+j, i] = Axis(fig[1+j, i], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} N_s")

            for k in 1:length(ϵϵ)
                label = L"ϵ=10^{%$(Int(log10(ϵϵ[k])))}"

                means = zeros(length(rr))
                lows = zeros(length(rr))
                highs = zeros(length(rr))

                for l in 1:length(rr)
                    skips = windows[l, k, i]
                    if length(skips) == 0 || j > length(skips) continue end
                    if length(skips[j]) == 0 continue end
                    σ = std(skips[j])
                    means[l] = max(mean(skips[j]), 1)   
                    low = min(max(means[l] - σ, 1), means[l])
                    lows[l] = low == 1 ? means[l] : low
                    highs[l] = max(max(means[l] + σ, 1), means[l])
                end

                real_data = findall(!=(0), means)

                if j == 1
                    lines!(axes[1, i], log10.(rr ./ rb), log10.(errors[:, k, i]), label=label)
                end
                lines!(axes[1+j, i], log10.(rr ./ rb)[real_data], log10.(means)[real_data])
                band!(axes[1+j, i], log10.(rr ./ rb)[real_data], log10.(lows)[real_data], log10.(highs)[real_data], alpha=0.3)
            end
        end
    end

    Legend(fig[max_K+2, :], axes[1, 1], orientation = :horizontal)

    for i in 1:length(q)
        Makie.ylims!(axes[1, i], -16, 0)
    end
    for i in 1:max_K
        for j in 1:length(q)
            hidexdecorations!(axes[i, j], grid = false)
        end
    end
    for i in 2:length(q)
        for j in 1:max_K+1
            hideydecorations!(axes[j, i], grid = false)
        end
    end

    for i in 1:max_K
        for j in 1:length(q)
            linkxaxes!(axes[i, j], axes[i+1, j])
        end
    end

    for i in 1:max_K+1
        for j in 1:length(q)-1
            linkyaxes!(axes[i, j], axes[i, j+1])
        end
    end

    for (i, title) in enumerate(titles)
        Label(fig[0, i], title, tellwidth=false, fontsize=20)
    end
    for i in 1:max_K
        Label(fig[i+1, 0], "Block $(i+1)", tellheight=false, fontsize=16, rotation = pi/2)
    end
    fig
end

function show_interpolation(ϵ, r, q)
    Δt = 3600.
    Nt = size(q)[2]

    α = 1e-6
    kg = 3.
    rb = 0.1
    constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

    positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
    windows = []

    Ib = make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
    I_int_bound = make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows)
    
    w1 = findfirst(>(0), I_int_bound[1,:])
    interpolation_limits = cumsum(vcat([w1], windows[1])) 
    interpolation_limits[end] = Nt

    plot_interpolation_limits = @. interpolation_limits * Δt ./ (8760 * 3600) 

    tt = @. (1:Nt) / (8760 * 3600) * Δt

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel=L"T \ (°C)", xlabel=L"t \ (years)",)
    lines!(ax, tt, Ib[1,:],          label = "Real")
    lines!(ax, tt, I_int_bound[1,:], label = "Interpolated", linewidth = 5, alpha = 0.6)

    vlines!(plot_interpolation_limits, linestyle = :dash, color = :green, label = "Interpolation \npoints")
    axislegend(position= :lt);
    fig
end

Nt = 8760*100
q1 = ones(2) * [sin(i/24) for i in 1:Nt]'
q2 = ones(2) * [sin(i/8760) for i in 1:Nt]'
q3 = ones(2, Nt)
q = [q1, q2, q3]
titles = [
    "High frequency",
    "Low frequency",
    "Constant"
]

high_osc_int = show_interpolation(1e-6, 10., q1[:, 1:5000])
low_osc_int = show_interpolation(5 * 1e-5, 10., q2[:, 1:20000])

save("interpolation_high.png", high_osc_int)
save("interpolation_low.png", low_osc_int)


fig_bound = generate_analysis_plot(q, use_bound=true, titles=titles)
fig_int = generate_analysis_plot(q, use_bound=false, titles=titles)

save("analysis_block_no_bound.png", fig_int)
save("analysis_block_bound.png", fig_bound)
