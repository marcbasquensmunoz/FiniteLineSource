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

function make_interpolated_simulation(;Nt, r, ϵ, positions, constants, q)
    block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);
    @unpack load_buffer, load_delays, qaux, expNin, expNout, ranges, expt, F, K_min, Kranges, HM = block

    Nb = size(K_min)[1]
    K = size(load_delays)[1]

    I_int = zeros(Nb, Nt)
    N_start = 1
    Ns = 1
    d2_count = 0
    current_q = zeros(Nb)

    q_der = zeros(3)
    d2Q = 0.
    dQ = 0.

    for nt in 1:Nt
        for i in 1:Nb
            current_q[i] = load_buffer[i][1]
        end
        @views for (q, b) in zip(q[:,nt], load_buffer)
            push!(b, q)
        end

        q_der[3] = q_der[2]
        q_der[2] = q_der[1]
        q_der[1] = current_q[1]

        d2Q = max(d2Q, abs((q_der[1] + q_der[3] - 2*q_der[2])/ Δt̃^2))
        dQ = max(dQ, abs((q_der[1] - q_der[2])/ Δt̃))

        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                push!(load_delays[i, j], qin)
                current_q[j] = qout
                @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
            end
            @views @. F[:, j] = expt * F[:, j] + qaux
        end


        if (nt-N_start)%Ns == 0 || nt == Nt
            for target in 1:Nb
                for source in 1:Nb
                    if source == target continue end
                    if K_min[source, target] == 0 || length(Kranges) < K_min[source, target] continue end
                    @inbounds range = Kranges[K_min[source, target]]
                    @inbounds @views I_int[target, nt] += dot(F[range, source], HM[range, target, source])
                end
                if nt == 1 continue end
                Nreal = nt == Nt ? Nt-N_start : Ns
                @views @. I_int[target, nt-Nreal+1:nt-1] += I_int[target, nt-Nreal] + (I_int[target, nt] - I_int[target, nt-Nreal]) * (1:Nreal-1) / Nreal
            end 
            d2_count += 1

            if d2_count >= 3 
                d2I = I_int[1, nt] + I_int[1, nt-2] - 2*I_int[1, nt-1]
                if d2I != 0 
                    Ns = Int(floor(sqrt(8ϵ/abs(d2I))))
                else 
                    Ns = 1
                end
                @show Ns
                        Q = q_der[1]
        N1 = block.N[1]
        B_disc = 1/(4*π^2*kg*r) * (Q * (2*log((N1+1)/(N2+1)) - log((N1*(N1+2))/(N2*(N2+2)))) / Δt̃^2 +  d2Q/2 * log(1+1/N1) + dQ/2 * ( 1/(N1 + N1^2)/Δt̃ + log(N1) - 2*log(1+N1) + log(2+N1)) + Q * 1/(2N1+3*N1^2+N1^3) / Δt̃ )
        @show 1/Δt̃ * sqrt(8*ϵ / abs(B_disc))
        d2Q = 0
        dQ = 0

                N_start = nt
                d2_count = 0
            else 
                Ns = 1
            end
        end
    end

    return I_int
end


function make_interpolated_simulation_with_bound(;Nt, r, ϵ, positions, constants, q, windows = [], K_int=1, use_bound=false)
    block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);
    @unpack load_buffer, load_delays, qaux, expNin, expNout, ranges, expt, F, K_min, Kranges, HM, ζ, N = block
    @unpack Δt̃, Δt, α, kg, rb = constants

    Nb = size(K_min)[1]
    K = size(load_delays)[1]

    ϵ = ϵ / K

    I_int = zeros(Nb, Nt)
    N_start = 1
    Ns = 0.
    current_q = zeros(Nb)

    q_in_der = (0., 0.)
    q_out_der = (0., 0.)
    dQ_in = 0.
    dQ_out = 0.
    Q_in = 0.
    Q_out = 0.
    Q_mean = 0.

    B_max = 0.

    N1 = length(N) >= K_int ? N[K_int] : return
    N2 = length(N) >= K_int+1 ? N[K_int+1] : N[K_int]

    source = 1
    target = 2


    C_1_in = (erf(r/rb / sqrt(4*N1*Δt̃)) - erf(r/rb / sqrt(4*(N1+1)*Δt̃))) / (4*π*kg*r * Δt̃) 
    C_1_out = -(erf(r/rb / sqrt(4*N2*Δt̃)) - erf(r/rb / sqrt(4*(N2+1)*Δt̃))) / (4*π*kg*r * Δt̃) 

    C_0_in = -( erf(r/rb / sqrt(4*N1*Δt̃)) + erf(r/rb / sqrt(4*(N1+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(N1+1)*Δt̃)) ) / (4*π*kg*r * Δt̃^2) 
    C_0_out = ( erf(r/rb / sqrt(4*N2*Δt̃)) + erf(r/rb / sqrt(4*(N2+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(N2+1)*Δt̃)) ) / (4*π*kg*r * Δt̃^2) 

    #derf = [(erf(r/rb / sqrt(4*Nb*Δt̃)) - erf(r/rb / sqrt(4*(Nb+1)*Δt̃))) / (4*π*kg*r * Δt̃)  for Nb in N]
    #d2erf = [( erf(r/rb / sqrt(4*Nb*Δt̃)) + erf(r/rb / sqrt(4*(Nb+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(Nb+1)*Δt̃)) ) / (4*π*kg*r * Δt̃^2) for Nb in N]

    F_1 = zeros(size(F))

    for nt in 1:Nt
        for i in 1:Nb
            current_q[i] = load_buffer[i][1]
        end
        @views for (q, b) in zip(q[:,nt], load_buffer)
            push!(b, q)
        end

        A = 0.
        B = 0.

        @views @. F_1 = F
        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                if i == K_int
                    A = qin
                    B = qout
                end
                push!(load_delays[i, j], qin)
                current_q[j] = qout
                @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
            end
            @views @. F[:, j] = expt * F[:, j] + qaux
        end
        
        if nt < N[1] || length(Kranges) == 0 continue end
        q_in_der = (A, q_in_der[1])
        q_out_der = (B, q_out_der[1])

        dQ_in = (q_in_der[1] - q_in_der[2]) / Δt
        dQ_out = (q_out_der[1] - q_out_der[2]) / Δt
        Q_in = q_in_der[1]
        Q_out = q_out_der[1]

        range = Kranges[K_min[source, target]]# ranges[K_int]#
        Q_eff = abs(Q_in) #max(abs(Q_in), abs(Q_out))
        osc_int = dot(F_1[range, source] .* (1 .- expt[range]) .^ 2, HM[range, target, source]) / Δt̃^2

        bounded_int = (log((N1+1)/(N2+1)) + 1/2 * log(N2*(N2+2)/(N1*(N1+2)))) / (2*π^2*kg*r*Δt̃^2) * Q_eff

        C_disc_1 = C_1_in * dQ_in + C_1_out * dQ_out 
        C_disc_2 = C_0_in * Q_in + C_0_out * Q_out 
        
        B_tot = abs(C_disc_1 + C_disc_2 + (use_bound ? bounded_int : osc_int))
        B_max = max(B_max, B_tot)
        B = 1/Δt̃ * sqrt(8ϵ / B_max)
        Ns = max(B == Inf ? 1 : Int(floor(B)), 1)

        if nt - N_start >= Ns || nt == Nt
            if nt > N[1] push!(windows, Ns) end
            for target in 1:Nb
                for source in 1:Nb
                    if source == target continue end
                    if K_min[source, target] == 0 || length(Kranges) < K_min[source, target] continue end
                    @inbounds range = Kranges[K_min[source, target]]
                    @inbounds @views I_int[target, nt] += dot(F[range, source], HM[range, target, source])
                end
                if nt == 1 continue end
                Nreal = nt - N_start 
                @views @. I_int[target, nt-Nreal+1:nt-1] += I_int[target, nt-Nreal] + (I_int[target, nt] - I_int[target, nt-Nreal]) * (1:Nreal-1) / Nreal
            end 
       
            N_start = nt
            Ns = 0
            B_max = 0.
        end
    end

    return I_int
end

ϵ = 1e-6
Δt = 3600.
Nt = 100

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

r = 1.
positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
#q = ones(2) * [0 + sin(i/24) for i in 1:Nt]'
q = ones(2, Nt)

windows = []

Ib = @time make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
#I_int = @time make_interpolated_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
I_int_bound = @time make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows, K_int=1, use_bound=false)

#=
maximum(abs.(Ib - I_int))
plot(log10.(abs.(Ib - I_int))[1,:])
=#
maximum(abs.(Ib - I_int_bound))

Plots.plot(log10.(abs.(Ib - I_int_bound))[1,:])

#=
plot(I_int[1,:])
plot!(Ib[1,:])
=#
Plots.plot(I_int_bound[1,:])
Plots.plot!(Ib[1,:])


Profile.clear()
@profile make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows, K_int=1, use_bound=false)
pprof()



function get_max_error(q, ϵ, r; block, use_bound)
    Δt = 3600.
    Nt = size(q)[2]

    α = 1e-6
    kg = 3.
    rb = 0.1
    constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

    positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
    windows = []

    Ib = make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
    I_int_bound = make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q, windows=windows, K_int=block, use_bound=use_bound)
    
    maximum(abs.(Ib - I_int_bound)), length(windows) == 0 ? 0. : mean(windows)
end

function generate_analysis_plot(q1, q2; block=1, use_bound=false)
    rr = 10. .^ (0:0.2:2)
    ϵϵ = 10. .^ (-2:-2:-10)

    err1 = zeros(length(rr), length(ϵϵ))
    W1 = zeros(length(rr), length(ϵϵ))

    err2 = zeros(length(rr), length(ϵϵ))
    W2 = zeros(length(rr), length(ϵϵ))

    for (i, r) in enumerate(rr)
        for (j, ϵ) in enumerate(ϵϵ)
            error1, Wm1 = get_max_error(q1, ϵ, r, block=block, use_bound=use_bound)
            error2, Wm2 = get_max_error(q2, ϵ, r, block=block, use_bound=use_bound)

            err1[i, j] = error1
            err2[i, j] = error2
            W1[i, j] = Wm1
            W2[i, j] = Wm2
        end
    end
    err1[err1 .< 1e-16] .= 1e-16
    err2[err2 .< 1e-16] .= 1e-16

    fig = Figure(size=(800, 600))

    ax11 = Axis(fig[1, 1], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} \ \Vert \epsilon  \Vert_{\infty}",)
    ax21 = Axis(fig[2, 1], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} N_s")
    ax12 = Axis(fig[1, 2], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} \ \Vert \epsilon  \Vert_{\infty}",)
    ax22 = Axis(fig[2, 2], xlabel=L"\log_{10} \tilde{r}", ylabel=L"\log_{10} N_s")

    for i in 1:length(ϵϵ)
        label = "ϵ=1e$(Int(log10(ϵϵ[i])))"
        lines!(ax11, log10.(rr ./ rb), log10.(err1[:, i]), label=label)
        lines!(ax21, log10.(rr ./ rb), log10.(W1[:, i]), label=label)
        lines!(ax12, log10.(rr ./ rb), log10.(err2[:, i]), label=label)
        lines!(ax22, log10.(rr ./ rb), log10.(W2[:, i]), label=label)
    end

    Legend(fig[1, 3], ax11)

    Makie.ylims!(ax11, -16, 0)
    Makie.ylims!(ax12, -16, 0)

    hidexdecorations!(ax11, grid = false)
    hidexdecorations!(ax12, grid = false)
    hideydecorations!(ax12, grid = false)
    hideydecorations!(ax22, grid = false)

    linkxaxes!(ax11, ax21)
    linkxaxes!(ax12, ax22)
    linkyaxes!(ax11, ax12)
    linkyaxes!(ax21, ax22)

    fig

    #save("analysis_block$(block)_bound$(use_bound).png", fig)
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
    interpolation_limits = cumsum(vcat([w1], windows)) 
    interpolation_limits[end] = Nt

    @. interpolation_limits *= Δt ./ (8760 * 3600) 

    tt = @. (1:Nt) / (8760 * 3600) * Δt

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel=L"T \ (°C)", xlabel=L"t \ (years)",)
    lines!(ax, tt, Ib[1,:],          label = "Real")
    lines!(ax, tt, I_int_bound[1,:], label = "Interpolated", linewidth = 3, alpha = 0.6)

    vlines!(interpolation_limits, linestyle = :dash, color = :green, label = "Interpolation \npoints")
    axislegend(position= :lt);
    fig
end

Nt = 8760*20
q1 = ones(2) * [0 + sin(i/24) for i in 1:Nt]'
q2 = ones(2, Nt)

show_interpolation(5 * 1e-4, 10., q1)
show_interpolation(1e-5, 10., q2)

generate_analysis_plot(q1, q2, block=1, use_bound=true)
generate_analysis_plot(q1, q2, block=1, use_bound=false)
generate_analysis_plot(q1, q2, block=2, use_bound=true)
generate_analysis_plot(q1, q2, block=2, use_bound=false)
