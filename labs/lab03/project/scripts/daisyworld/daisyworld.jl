using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, CairoMakie

include(srcdir("daisyworld.jl"))

model = daisyworld()
daisycolor(a::Daisy) = a.breed

plotkwargs = (
    agent_color=daisycolor, agent_size=20, agent_marker='✿',
    heatarray=:temperature, heatkwargs=(colorrange=(-20, 60),)
)

plt1, _ = abmplot(model; plotkwargs...)
save(plotsdir("daisy_step001.png"), plt1)

step!(model, 5)
plt2, _ = abmplot(model; heatarray=model.temperature, plotkwargs...)
save(plotsdir("daisy_step005.png"), plt2)

step!(model, 40)
plt3, _ = abmplot(model; heatarray=model.temperature, plotkwargs...)
save(plotsdir("daisy_step045.png"), plt3)
