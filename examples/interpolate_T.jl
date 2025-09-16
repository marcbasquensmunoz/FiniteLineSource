using FiniteLineSource
using FiniteLineSource: PointSource
using Parameters
using LinearAlgebra

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


function make_interpolated_simulation_with_bound(;Nt, r, ϵ, positions, constants, q)
    block = prepare_containers(PointToPoint(r = r), positions, ϵ, Nt, constants, compute_first_block=false);
    @unpack load_buffer, load_delays, qaux, expNin, expNout, ranges, expt, F, K_min, Kranges, HM, ζ, N = block

    Nb = size(K_min)[1]
    K = size(load_delays)[1]

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

    B_max = 0.

    N1 = N[1]
    N2 = length(N) >= 2 ? N[2] : N[1]

    for nt in 1:Nt
        for i in 1:Nb
            current_q[i] = load_buffer[i][1]
        end
        @views for (q, b) in zip(q[:,nt], load_buffer)
            push!(b, q)
        end

        A = 0.
        B = 0.

        for j in 1:Nb
            qaux .= 0.
            for i in 1:K
                qin = current_q[j]
                qout = i == K ? 0. : load_delays[i, j][1]
                A = qin
                B = qout
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

        range = Kranges[K_min[1, 2]]
        loc = dot(F[range, 1] .* ζ[range] .^4, HM[range, 2, 1]) / Δt̃^2 * max(abs(Q_in), abs(Q_out))
        comp = (log(N1+1) - log(N2+1) + 1/2 * log(N2*(N2+2)) - 1/2*log(N1*(N1+2))) / (2*π^2*kg*r*Δt̃^2) * max(abs(Q_in), abs(Q_out))
        #@show loc, comp

        #B_disc_1 = 1 / (4π*kg*r) * ( dQ_in * erf(r/rb / sqrt(4*N1*Δt̃)) - dQ_out * erf(r/rb / sqrt(4*N2*Δt̃)) )
        #B_disc_2 = 1 / ((4π)^(3/2)*kg) * (Q_in * exp(-(r/rb)^2 / (4*(N1*Δt̃))) / (N1*Δt̃)^(3/2) - Q_out * exp(-(r/rb)^2 / (4*(N2*Δt̃))) / (N2*Δt̃)^(3/2) )
        #B_tot = abs(loc + B_disc_1 - B_disc_2)

        C_disc_1 = 1 / (4π*kg*r) * ( dQ_in * (erf(r/rb / sqrt(4*N1*Δt̃)) - erf(r/rb / sqrt(4*(N1+1)*Δt̃)))/Δt̃ -  dQ_out * (erf(r/rb / sqrt(4*N2*Δt̃)) - erf(r/rb / sqrt(4*(N2+1)*Δt̃)))/Δt̃ )
        C_disc_2 = - 1 / (4π*kg*r) * ( 
              Q_in * ( erf(r/rb / sqrt(4*N1*Δt̃)) + erf(r/rb / sqrt(4*(N1+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(N1+1)*Δt̃)) ) / Δt̃^2 
            - Q_out * ( erf(r/rb / sqrt(4*N2*Δt̃)) + erf(r/rb / sqrt(4*(N2+2)*Δt̃)) - 2*erf(r/rb / sqrt(4*(N2+1)*Δt̃)) ) / Δt̃^2 
        )

        B_tot = abs(C_disc_1 + C_disc_2 + loc)
        #B_tot = abs(C_disc_1 + C_disc_2 + comp)

        B_max = max(B_max, B_tot)
        B = 1/Δt̃ * sqrt(8ϵ / B_max)
        Ns = max(B == Inf ? 1 : Int(floor(B)), 1)

        #@show Ns, 1/Δt̃ * sqrt(8ϵ / abs(C))
        #@show B_disc_1, B_disc_2, loc
        
        if nt - N_start >= Ns || nt == Nt
            @show nt, Ns
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

ϵ = 1e-8
Δt = 3600.
Nt = 8760

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

r = 1.
positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
q = ones(2) * [sin(i/24) for i in 1:Nt]'

Ib = @time make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
#I_int = @time make_interpolated_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
I_int_bound = @time make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)

#=
maximum(abs.(Ib - I_int))
plot(log10.(abs.(Ib - I_int))[1,:])
=#
maximum(abs.(Ib - I_int_bound))
plot(log10.(abs.(Ib - I_int_bound))[1,:])

#=
plot(I_int[1,:])
plot!(Ib[1,:])
=#
plot(I_int_bound[1,:])
plot!(Ib[1,:])


function get_max_error(ϵ, r)
    Δt = 3600.
    Nt = 8760

    α = 1e-6
    kg = 3.
    rb = 0.1
    constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

    positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
    q = ones(2) * [sin(i/24) for i in 1:Nt]'
    #q = ones(2, Nt)

    Ib = @time make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
    I_int_bound = @time make_interpolated_simulation_with_bound(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)

    maximum(abs.(Ib - I_int_bound))
end

using CairoMakie

rr = [1., 2., 5., 10., 20., 50., 100.]
ϵϵ = 10. .^ (-2:-2:-10)
err = [get_max_error(ϵ, r) for r in rr, ϵ in ϵϵ]
err[err .< 1e-16] .= 1e-16

fig = Figure()
ax = Axis(fig[1, 1])
for i in 1:length(ϵϵ)
    lines!(ax, log10.(rr), log10.(err[:, i]), label = "ϵ=1e$(Int(log10(ϵϵ[i])))")
end

#Makie.xlims!(ax, 1, 10)
Makie.ylims!(ax, -16, 0)
axislegend(position= :rt);

fig
