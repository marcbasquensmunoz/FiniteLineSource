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

    for nt in 1:Nt
        for i in 1:Nb
            current_q[i] = load_buffer[i][1]
        end
        @views for (q, b) in zip(q[:,nt], load_buffer)
            push!(b, q)
        end

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
                N_start = nt
                d2_count = 0
            else 
                Ns = 1
            end
        end
    end

    return I_int
end


ϵ = 1e-6
Δt = 3600.
Nt = 8760*20

α = 1e-6
kg = 3.
rb = 0.1
constants = Constants(Δt=Δt, α=α, kg=kg, rb=rb)

r = 10.
positions = [PointSource(0., 0., 0., rb), PointSource(r, 0., 0., rb)]
q = ones(2) * [sin(i/30) for i in 1:Nt]'

Ib = @btime make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
I_int = @btime make_interpolated_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)

maximum(abs.(Ib - I_int))
plot(log10.(abs.(Ib - I_int))[1,:])

plot(I_int[1,:])
plot!(Ib[1,:])

Profile.Allocs.clear()
@time Profile.Allocs.@profile sample_rate=1 make_full_simulation(Nt=Nt, r=r, ϵ=ϵ, positions=positions, constants=constants, q=q)
PProf.Allocs.pprof(from_c=false)


#=
Qp = maximum(abs.(diff(q[1,:]))) / Δt̃
Q = 1
N1 = 886
N2 = 8760

B_disc = 1/(4*π^2*kg*r * Δt̃^2) * ( 
    Q * (2*log((N1+1)/(N2+1)) - log((N1*(N1+2))/(N2*(N2+2))))
    + (Δt̃ * Qp * log(1+1/N1) - Q * log((N1+1)^2/(N1*(N1+2))))
    - (Δt̃ * Qp * log(1+1/N2) - Q * log((N2+1)^2/(N2*(N2+2))))
) 
1/Δt̃ * sqrt(8*ϵ / B_disc)


B_disc = 1/(2*π^2*kg*r * Δt̃^2) * ( 
    -Q * (log(N1+2) - 2*log(N1+1) + log(N1))
    + Qp *Δt̃ * (log(1+1/N1) + log(1+1/N2))
) 
1/Δt̃ * sqrt(8*ϵ / B_disc)



B_disc = 1/(4*π^2*kg*r * Δt̃^2) * ( 
    2 * Q * (2*log((N1+1)/(N2+1)) - log((N1*(N1+2))/(N2*(N2+2))))
    + Qp * Δt̃ / 2 * log((N2*(N1+1))/(N1*(N2+1)))
) 
1/Δt̃ * sqrt(8*ϵ / B_disc)
=#
