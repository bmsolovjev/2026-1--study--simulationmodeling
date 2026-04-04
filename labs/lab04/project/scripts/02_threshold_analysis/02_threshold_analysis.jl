using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, Statistics

include(srcdir("sir_model.jl"))

function find_peak(beta; seed=42)
    model = initialize_sir(;
        Ns=[1000, 1000, 1000],
        β_und=fill(beta, 3),
        β_det=fill(beta/10, 3),
        infection_period=14,
        detection_time=7,
        death_rate=0.02,
        reinfection_probability=0.1,
        Is=[0, 0, 1],
        seed=seed,
        n_steps=100,
    )

    infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
    peak = 0.0

    for step = 1:100
        Agents.step!(model, 1)
        frac = infected_frac(model)
        if frac > peak
            peak = frac
        end
    end

    return peak
end

beta_range = 0.05:0.01:0.7
peaks = Float64[]

for β in beta_range
    push!(peaks, find_peak(β))
end

df = DataFrame(beta=beta_range, peak=peaks)

threshold_beta = minimum(df.beta[df.peak .> 0.05])

println("Минимальное β для эпидемии (пик > 5%): ", round(threshold_beta, digits=3))

γ = 1 / 14
R₀_threshold = 1
β_theoretical = γ * R₀_threshold
println("Теоретический порог (R₀ = 1): β = ", round(β_theoretical, digits=3))

plot(df.beta, df.peak .* 100, linewidth=2, label="Пик инфицированных")
vline!([threshold_beta], linestyle=:dash, label="Практический порог (5%)", color=:red)
vline!([β_theoretical], linestyle=:dot, label="Теоретический порог (R₀=1)", color=:blue)

xlabel!("Коэффициент заразности β")
ylabel!("Пиковая доля инфицированных (%)")
title!("Порог возникновения эпидемии")
savefig(plotsdir("sir_threshold_analysis.png"))

println("Результаты сохранены в plots/sir_threshold_analysis.png")
