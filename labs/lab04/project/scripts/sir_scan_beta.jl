# # Сканирование параметра β в модели SIR
# Исследуем влияние коэффициента заразности на динамику эпидемии.
using DrWatson
@quickactivate "project"

using Agents, DataFrames, Plots, CSV, Random

# Подключаем файл с определением модели.
include(srcdir("sir_model.jl"))

# ## Функция для запуска одного эксперимента
# Принимает параметры и возвращает ключевые метрики.
function run_experiment(p)
    # Создаём β_und и β_det на основе скалярного beta
    beta = p[:beta]
    β_und = fill(beta, 3)
    β_det = fill(beta/10, 3)

    # Инициализируем модель
    model = initialize_sir(;
        Ns = p[:Ns],
        β_und = β_und,
        β_det = β_det,
        infection_period = p[:infection_period],
        detection_time = p[:detection_time],
        death_rate = p[:death_rate],
        reinfection_probability = p[:reinfection_probability],
        Is = p[:Is],
        seed = p[:seed],
        n_steps = p[:n_steps],
    )

    # Вспомогательная функция для доли инфицированных
    infected_fraction(model) =
        count(a.status == :I for a in allagents(model)) / nagents(model)

    peak_infected = 0.0

    # Ручной шаг симуляции (безопасный)
    for step = 1:p[:n_steps]
        agent_ids = collect(allids(model))
        for id in agent_ids
            agent = try
                model[id]
            catch
                nothing
            end
            if agent !== nothing
                sir_agent_step!(agent, model)
            end
        end

        frac = infected_fraction(model)
        if frac > peak_infected
            peak_infected = frac
        end
    end

    final_infected = infected_fraction(model)
    final_recovered = count(a.status == :R for a in allagents(model)) / nagents(model)
    total_deaths = sum(p[:Ns]) - nagents(model)

    return (
        peak = peak_infected,
        final_inf = final_infected,
        final_rec = final_recovered,
        deaths = total_deaths,
    )
end

# ## Параметры эксперимента
# Диапазон значений β и повторы для разных seed.
beta_range = 0.1:0.1:1.0
seeds = [42, 43, 44]

# ## Генерация списка параметров
# Создаём все комбинации β и seed.
params_list = []
for b in beta_range
    for s in seeds
        push!(
            params_list,
            Dict(
                :beta => b,
                :Ns => [1000, 1000, 1000],
                :infection_period => 14,
                :detection_time => 7,
                :death_rate => 0.02,
                :reinfection_probability => 0.1,
                :Is => [0, 0, 1],
                :seed => s,
                :n_steps => 100,
            ),
        )
    end
end

# ## Запуск экспериментов
# Выполняем все прогоны и собираем результаты.
results = []
for params in params_list
    data = run_experiment(params)
    push!(results, merge(params, Dict(pairs(data))))
    println("Завершён эксперимент с beta = $(params[:beta]), seed = $(params[:seed])")
end

# ## Сохранение результатов
# Сохраняем данные всех прогонов в CSV.
df = DataFrame(results)
CSV.write(datadir("beta_scan_all.csv"), df)

# ## Усреднение по повторам
# Усредняем метрики для каждого значения β.
using Statistics
grouped = combine(
    groupby(df, [:beta]),
    :peak => mean => :mean_peak,
    :final_inf => mean => :mean_final_inf,
    :deaths => mean => :mean_deaths,
)

# ## Визуализация результатов
# Строим график зависимости метрик от β.
plot(
    grouped.beta,
    grouped.mean_peak,
    label = "Пик эпидемии",
    xlabel = "Коэффициент заразности β",
    ylabel = "Доля инфицированных",
    marker = :circle,
    linewidth = 2,
)
plot!(
    grouped.beta,
    grouped.mean_final_inf,
    label = "Конечная доля инфицированных",
    marker = :square,
)
plot!(grouped.beta, grouped.mean_deaths ./ 3000, label = "Доля умерших", marker = :diamond)

# Сохраняем график.
savefig(plotsdir("beta_scan.png"))

println("Результаты сохранены в data/beta_scan_all.csv и plots/beta_scan.png")