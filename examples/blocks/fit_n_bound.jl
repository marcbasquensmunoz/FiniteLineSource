using Flux, Statistics, ProgressMeter, Random
using JLD2
using Plots, Parameters


vσ = [1., 5., 10., 20., 50., 100., 200., 400., 700., 1000.]
vb = [2., 4., 6., 10., 20., 35., 50., 75., 100., 200., 500., 1000.]
va = [1., 3., 5., 8., 10., 15., 25., 35., 50., 75., 100.]
vϵ = [10. ^k for k in -4:-2:-12]

N = length(vσ)*length(vb)*length(va)*length(vϵ)

x_train = zeros(Float32, 4, N)
y_train = zeros(Float32, N)

index = 1
for (i, ϵ) in enumerate(vϵ)
    for (j, σ) in enumerate(vσ)
        for (k, A) in enumerate(va)
            for (l, B) in enumerate(vb)
                b = σ * B
                if k == 1 || σ + 2 > σ * A
                    a = min(σ + max(2, (b-σ) * 0.01), b)
                else 
                    a = σ * A
                end
                x_train[:, index] = [-log10(ϵ), log(σ), log(a), log(b)]
                #y_train[index] = min_n[i, j, k, l]
                index += 1
            end
        end
    end
end

mins = zeros(4)
norms = zeros(4)

for i in 1:4
    mins[i] = minimum(x_train[i, :])
    norms[i] = maximum(x_train[i, :]) - mins[i]
    @. x_train[i, :] = (x_train[i, :] - mins[i]) / norms[i]
end

data = [(Float32.(x_train[:, i]), y_train[i]) for i in 1:N]

k = 32
model = Chain(
    Dense(4 => k, relu),   
    Dense(k => k, relu),   
    Dense(k => 1)
)

loss(model, x, y) = mean(abs2.(model(x) .- y))
opt_state = Flux.setup(Adam(0.001), model)

for epoch in 1:20
    if epoch % 10 == 0
        @info "Epoch $epoch"
    end
    Flux.train!(loss, model, data, opt_state)  
end

eval(ϵ, σ, a, b, model, mins, norms) = model(Float32[(-log10(ϵ) - mins[1]) / norms[1], (log(σ) - mins[2]) / norms[2], (log(a) - mins[3]) / norms[3], (log(b) - mins[4]) / norms[4]])[1]
#eval(ϵ, σ, a, b, model) = model(Float32[-log10(ϵ), log(σ), log(a), log(b)])[1]

bb = 10.:1.:1000.
res = eval.(Ref(nmodel), 1e-6, 1., 8., bb) 
#res = eval.(1e-6, 1., 3., bb, Ref(model), Ref(mins), Ref(norms))
plot(bb, res)   
scatter!(vb, min_n[2, 1, 4, :])


# Save the model
#=
model_state = Flux.state(model)
jldsave("N_bound.jld2"; model_state)
=#
#= Load the model
model_state = JLD2.load("N_bound.jld2", "model_state")
Flux.loadmodel!(model, model_state)
=#