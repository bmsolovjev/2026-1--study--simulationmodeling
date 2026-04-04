# # Многокритериальная оптимизация с ограничениями (увеличенная нижняя граница β)
using DrWatson
@quickactivate "project"

using BlackBoxOptim, Random, Statistics
using DrWatson: @save

include(srcdir("sir_model.jl"))

# ## Целевая функция с ограничениями
function cost_with_constraint(x)
    β_und = x[1]
    detection_time = round(Int, x[2])
    death_rate = x[3]

    replicates = 5
    peak_vals = Float64[]
    dead_vals = Float64[]

    for rep in 1:replicates
        model = initialize_sir(;
            Ns=[1000, 1000, 1000],
            β_und=fill(β_und, 3),
            β_det=fill(β_und/10, 3),
            infection_period=14,
            detection_time=detection_time,
            death_rate=death_rate,
            reinfection_probability=0.1,
            Is=[0, 0, 1],
            seed=42 + rep,
            n_steps=100,
        )

        infected_frac(model) = count(a.status == :I for a in allagents(model)) / nagents(model)
        peak = 0.0

        for step in 1:100
            Agents.step!(model, 1)
            frac = infected_frac(model)
            if frac > peak
                peak = frac
            end
        end

        deaths = 3000 - nagents(model)

        push!(peak_vals, peak)
        push!(dead_vals, deaths / 3000)
    end

    mean_peak = mean(peak_vals)
    mean_deaths = mean(dead_vals)

    # Штраф за превышение пика > 30%
    if mean_peak > 0.3
        return (Inf, mean_deaths)
    end

    return (mean_peak, mean_deaths)
end

# ## Запуск оптимизации с увеличенной нижней границей β
result = bboptimize(
    cost_with_constraint,
    Method = :borg_moea,
    FitnessScheme = ParetoFitnessScheme{2}(is_minimizing=true),
    SearchRange = [
        (0.3, 1.0),    # β_und — увеличенная нижняя граница с 0.1 до 0.3
        (3.0, 14.0),   # detection_time
        (0.01, 0.1),   # death_rate
    ],
    NumDimensions = 3,
    MaxTime = 120,
    TraceMode = :compact,
)

best = best_candidate(result)
fitness = best_fitness(result)

println("Оптимальные параметры (при пике < 30%):")
println("β_und = ", round(best[1], digits=3))
println("Время выявления = ", round(Int, best[2]), " дней")
println("Смертность = ", round(best[3], digits=4))
println()
println("Достигнутые показатели:")
println("Пик заболеваемости: ", round(fitness[1] * 100, digits=2), "%")
println("Доля умерших: ", round(fitness[2] * 100, digits=4), "%")

result_dict = Dict("best" => best, "fitness" => fitness)
save(datadir("constrained_optimization_result.jld2"), result_dict)