using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots

include(srcdir("sir_model.jl"))

function quarantine_model(; quarantine_threshold=0.1, seed=42)
    C = 3
    migration_rates = fill(0.1 / (C-1), C, C)
    for i in 1:C
        migration_rates[i, i] = 0.9
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

    infected_frac(model, city) = count(a.status == :I for a in allagents(model) if a.pos == city) / model.Ns[city]
    quarantined = false

    times = Int[]
    I_total = Float64[]
    I_cities = [[], [], []]

    for step in 1:150
        Agents.step!(model, 1)

        push!(times, step)
        push!(I_total, count(a.status == :I for a in allagents(model)) / nagents(model))

        for city in 1:3
            push!(I_cities[city], infected_frac(model, city))
        end

        if !quarantined
            for city in 1:3
                if infected_frac(model, city) > quarantine_threshold

                    model.migration_rates[city, :] .= 0
                    model.migration_rates[city, city] = 1
                    quarantined = true
                    println("Город $city закрыт на карантин на шаге $step")
                    break
                end
            end
        end
    end

    return (times=times, I_total=I_total, I_cities=I_cities)
end

function no_quarantine_model(seed=42)
    C = 3
    migration_rates = fill(0.1 / (C-1), C, C)
    for i in 1:C
        migration_rates[i, i] = 0.9
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

    times = Int[]
    I_total = Float64[]

    for step in 1:150
        Agents.step!(model, 1)
        push!(times, step)
        push!(I_total, count(a.status == :I for a in allagents(model)) / nagents(model))
    end

    return (times=times, I_total=I_total)
end

result_q = quarantine_model(quarantine_threshold=0.1)
result_nq = no_quarantine_model()

plot(result_nq.times, result_nq.I_total .* 100,
     label="Без карантина", linewidth=2,
     xlabel="Дни", ylabel="Доля инфицированных (%)",
     title="Эффективность карантинных мер")

plot!(result_q.times, result_q.I_total .* 100,
      label="С карантином (порог 10%)", linewidth=2)

savefig(plotsdir("sir_quarantine_effect.png"))

println("График сохранён в plots/sir_quarantine_effect.png")
