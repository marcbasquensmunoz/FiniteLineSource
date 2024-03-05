using FiniteLineSource
using Plots 

params = LevelSetParams(D1=8., H1=10., D2=0., H2=20., σ=0., rb=1.)

f, r = level_set_length(params)

plot(r, f.(r))