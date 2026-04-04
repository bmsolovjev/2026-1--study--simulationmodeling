# # Исследование влияния миграции на скорость распространения
using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, Statistics

include(srcdir("sir_model.jl"))

# ## Функция для анализа времени пика
function analyze_migration(intensity; seed=42)
    C = 3
    migration_rates = fill(intensity / (C-1), C, C)
    for i in 1:C
        migration_rates[i, i] = 1 - intensity
    end

    model = initialize_sir(;
        Ns=[1000, 1000, 1000],
        β_und=fill(0.5, 3),
        β_det=fill(0.05, 3),
        infection_period=14,
        detection_time=7,
        death_rate=0.02,
        reinfection_probability=0.1,
        Is=[1, 0, 0],
        seed=seed,
        migration_rates=migration_rates,
        n_steps=150,
    )

    infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
    peak = 0.0
    peak_time = 0

    for step in 1:150
        Agents.step!(model, 1)
        frac = infected_frac(model)
        if frac > peak
            peak = frac
            peak_time = step
        end
    end

    return (peak_time=peak_time, peak_value=peak)
end

# ## Сканирование интенсивности миграции
intensities = 0.0:0.05:0.5
seeds = [42, 43, 44]

results = []
for inten in intensities
    for s in seeds
        data = analyze_migration(inten; seed=s)
        push!(results, (intensity=inten, seed=s, peak_time=data.peak_time, peak_value=data.peak_value))
    end
end

df = DataFrame(results)
grouped = combine(groupby(df, :intensity),
                  :peak_time => mean => :mean_peak_time,
                  :peak_value => mean => :mean_peak_value)

# ## Визуализация
p1 = plot(grouped.intensity, grouped.mean_peak_time,
          marker=:circle, linewidth=2,
          xlabel="Интенсивность миграции",
          ylabel="Время до пика (дни)",
          title="Влияние миграции на скорость распространения")

p2 = plot(grouped.intensity, grouped.mean_peak_value .* 100,
          marker=:square, linewidth=2,
          xlabel="Интенсивность миграции",
          ylabel="Пиковая заболеваемость (%)",
          title="Влияние миграции на масштаб эпидемии")

plot(p1, p2, layout=(2,1), size=(600, 600))
savefig(plotsdir("sir_migration_analysis.png"))

min_time_idx = argmin(grouped.mean_peak_time)
println("Минимальное время до пика: ", round(grouped.mean_peak_time[min_time_idx], digits=0), " дней")
println("При интенсивности миграции: ", grouped.intensity[min_time_idx])