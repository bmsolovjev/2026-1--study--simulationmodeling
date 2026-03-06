using DrWatson
quickactivate(dirname(@__DIR__), "project")

using DifferentialEquations
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots
using Statistics
using FFTW
using LinearAlgebra
using CSV

script_name = splitext(basename(PROGRAM_FILE))[1]
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))


# ## МОДЕЛЬ ЛОТКИ-ВОЛЬТЕРРЫ (ХИЩНИК-ЖЕРТВА)

"""
Модель Лотки-Вольтерры (хищник-жертва)

Данная модель описывает динамику взаимодействия двух видов:
  - Жертвы (x) — популяция, служащая пищей для хищников
  - Хищники (y) — популяция, питающаяся жертвами

Система дифференциальных уравнений:
  dx/dt = α·x - β·x·y    (1)  — рост жертв ограничен хищниками
  dy/dt = δ·x·y - γ·y    (2)  — рост хищников зависит от наличия жертв

Физический смысл параметров:
  α [время⁻¹] — коэффициент естественного прироста жертв
                (рождаемость минус естественная смертность)
  β [особь⁻¹·время⁻¹] — коэффициент эффективности охоты хищников
                (вероятность встречи и поедания)
  δ [особь⁻¹·время⁻¹] — коэффициент конверсии биомассы
                (сколько хищников может прокормить одна жертва)
  γ [время⁻¹] — коэффициент естественной смертности хищников

Историческая справка:
  Модель была независимо предложена Альфредом Лоткой (1925) и
  Вито Вольтеррой (1926). Вольтерра использовал её для объяснения
  колебаний уловов рыб в Адриатическом море после Первой мировой войны.
"""
function lotka_volterra!(du, u, p, t)
    x, y = u               # x - жертвы, y - хищники
    α, β, δ, γ = p          # параметры модели

    @inbounds begin
        du[1] = α*x - β*x*y      # dx/dt: изменение популяции жертв
        du[2] = δ*x*y - γ*y      # dy/dt: изменение популяции хищников
    end
    nothing
end


# ## АНАЛИТИЧЕСКИЕ СВОЙСТВА МОДЕЛИ

"""
Анализ стационарных точек модели Лотки-Вольтерры

Система имеет две стационарные точки:

1. Тривиальная точка: (x, y) = (0, 0)
   — соответствует вымиранию обоих видов.
   Матрица Якоби: J = [α  0; 0 -γ]
   Собственные числа: λ₁ = α (>0), λ₂ = -γ (<0)
   → седловая точка (неустойчива)

2. Нетривиальная точка: (x*, y*) = (γ/δ, α/β)
   — точка сосуществования видов.
   Матрица Якоби: J = [0  -βx*; δy*  0]
   Собственные числа: чисто мнимые
   → центр (устойчивые колебания)

В окрестности нетривиальной стационарной точки решение имеет вид:
  x(t) ≈ x* + A·cos(ωt + φ)
  y(t) ≈ y* + B·sin(ωt + φ)

Где частота колебаний: ω = √(αγ)
Период колебаний: T = 2π/√(αγ)
"""
function analyze_fixed_points(p)
    α, β, δ, γ = p

    # Стационарная точка сосуществования
    x_star = γ / δ
    y_star = α / β

    # Частота и период колебаний
    ω = sqrt(α * γ)
    T = 2π / ω

    # Проверка устойчивости (линеаризация)
    J = [0  -β*x_star; δ*y_star  0]
    eigenvalues = eigvals(J)

    return (;
        x_star = x_star,
        y_star = y_star,
        frequency = ω,
        period = T,
        eigenvalues = eigenvalues,
        is_stable = all(real(eigenvalues) .<= 0)
    )
end


# ## НАБОРЫ ПАРАМЕТРОВ ДЛЯ ВЫЧИСЛЕНИЙ

"""
Наборы параметров модели Лотки-Вольтерры

Каждый набор соответствует различным экологическим сценариям:
  1. Классический сценарий — базовые параметры из учебников
  2. Быстрый сценарий — высокая скорость размножения жертв
  3. Медленный сценарий — низкая скорость размножения жертв
"""
parameter_sets = [
    # [α, β, δ, γ]
    [0.1, 0.02, 0.01, 0.3],    # Набор 1: Классический (умеренные колебания)
    [0.2, 0.02, 0.01, 0.3],    # Набор 2: Быстрый рост жертв (высокая частота)
    [0.05, 0.02, 0.01, 0.3],   # Набор 3: Медленный рост жертв (низкая частота)
    [0.1, 0.04, 0.02, 0.3],    # Набор 4: Высокая эффективность хищников
    [0.1, 0.01, 0.005, 0.15]   # Набор 5: Сбалансированный (устойчивые колебания)
]

parameter_names = [
    "Классический (умеренные колебания)",
    "Быстрый рост жертв (высокая частота)",
    "Медленный рост жертв (низкая частота)",
    "Высокая эффективность хищников",
    "Сбалансированный"
]

# Начальные условия (одинаковые для всех экспериментов)
u0_lv = [40.0, 9.0]  # [жертвы, хищники]

# Временные параметры
tspan_lv = (0.0, 200.0)   # длительность симуляции
dt_lv = 0.01              # шаг интегрирования


# ## ФУНКЦИЯ ДЛЯ РЕШЕНИЯ МОДЕЛИ С ЗАДАННЫМИ ПАРАМЕТРАМИ

"""
    solve_lotka_volterra(p, u0, tspan; kwargs...)

Решает систему Лотки-Вольтерры с заданными параметрами.

Аргументы:
  p — вектор параметров [α, β, δ, γ]
  u0 — начальные условия [x0, y0]
  tspan — интервал времени (t_start, t_end)

Возвращает структуру с решением и проанализированными данными.
"""
function solve_lotka_volterra(p, u0, tspan; saveat=0.1, reltol=1e-8, abstol=1e-10)
    # Создание и решение задачи
    prob = ODEProblem(lotka_volterra!, u0, tspan, p)
    sol = solve(prob,
                Tsit5(),
                dt = dt_lv,
                reltol = reltol,
                abstol = abstol,
                saveat = saveat,
                dense = true)

    # Подготовка данных
    df = DataFrame()
    df[!, :t] = sol.t
    df[!, :prey] = [u[1] for u in sol.u]
    df[!, :predator] = [u[2] for u in sol.u]

    # Расчет производных
    α, β, δ, γ = p
    df[!, :dprey_dt] = α .* df.prey .- β .* df.prey .* df.predator
    df[!, :dpredator_dt] = δ .* df.prey .* df.predator .- γ .* df.predator

    # Относительные изменения
    df[!, :prey_rel_change] = df.dprey_dt ./ df.prey .* 100
    df[!, :predator_rel_change] = df.dpredator_dt ./ df.predator .* 100

    return (sol = sol, df = df)
end


# ## АНАЛИЗ РЕЗУЛЬТАТОВ ДЛЯ ВСЕХ НАБОРОВ ПАРАМЕТРОВ

println("\n" * "="^80)
println(" " * "МОДЕЛЬ ЛОТКИ-ВОЛЬТЕРРЫ: АНАЛИЗ НАБОРОВ ПАРАМЕТРОВ")
println("="^80)

# Хранение результатов всех экспериментов
results = []

for (idx, p) in enumerate(parameter_sets)
    α, β, δ, γ = p

    println("\n" * "-"^60)
    println("ЭКСПЕРИМЕНТ #$idx: $(parameter_names[idx])")
    println("-"^60)

    # Аналитический анализ стационарных точек
    fp_analysis = analyze_fixed_points(p)
    println("\n📊 Аналитический анализ:")
    println("  Стационарная точка: (x*, y*) = (",
            round(fp_analysis.x_star, digits=2), ", ",
            round(fp_analysis.y_star, digits=2), ")")
    println("  Теоретическая частота колебаний: ω = ", round(fp_analysis.frequency, digits=4))
    println("  Теоретический период: T = ", round(fp_analysis.period, digits=2))

    # Численное решение
    println("\n🔬 Численное моделирование:")

    @time begin
        result = solve_lotka_volterra(p, u0_lv, tspan_lv)
        sol, df = result.sol, result.df
    end

    # Статистический анализ
    prey_min = minimum(df.prey)
    prey_max = maximum(df.prey)
    prey_mean = mean(df.prey)
    predator_min = minimum(df.predator)
    predator_max = maximum(df.predator)
    predator_mean = mean(df.predator)

    println("\n📈 Статистика популяций:")
    println("  Жертвы: min = ", round(prey_min, digits=2),
            ", max = ", round(prey_max, digits=2),
            ", среднее = ", round(prey_mean, digits=2))
    println("  Хищники: min = ", round(predator_min, digits=2),
            ", max = ", round(predator_max, digits=2),
            ", среднее = ", round(predator_mean, digits=2))

    # Спектральный анализ (нахождение доминирующей частоты)
    function find_dominant_frequency(signal, dt)
        n = length(signal)
        spectrum = abs.(rfft(signal .- mean(signal)))
        freq = rfftfreq(n, 1/dt)

        # Ищем максимум, исключая нулевую частоту
        idx_max = argmax(spectrum[2:end]) + 1
        return freq[idx_max], spectrum[idx_max]
    end

    dom_freq_prey, _ = find_dominant_frequency(df.prey, dt_lv)
    dom_freq_pred, _ = find_dominant_frequency(df.predator, dt_lv)

    println("\n🎵 Спектральный анализ:")
    println("  Доминирующая частота (жертвы): ", round(dom_freq_prey, digits=4))
    println("  Экспериментальный период (жертвы): ", round(1/dom_freq_prey, digits=2))
    println("  Доминирующая частота (хищники): ", round(dom_freq_pred, digits=4))
    println("  Экспериментальный период (хищники): ", round(1/dom_freq_pred, digits=2))

    # Проверка соответствия теории
    period_error = abs(1/dom_freq_prey - fp_analysis.period) / fp_analysis.period * 100
    println("\n🎯 Соответствие теории:")
    println("  Отклонение периода от теоретического: ", round(period_error, digits=2), "%")

    # Сохранение результатов
    push!(results, (;
        params = p,
        name = parameter_names[idx],
        df = df,
        sol = sol,
        stats = (min_prey=prey_min, max_prey=prey_max, mean_prey=prey_mean,
                 min_pred=predator_min, max_pred=predator_max, mean_pred=predator_mean),
        fp = fp_analysis,
        dom_freq = (prey=dom_freq_prey, predator=dom_freq_pred)
    ))

    # Построение графика для текущего набора параметров
    plt = plot(df.t, [df.prey df.predator],
               label=[L"Жертвы (x)" L"Хищники (y)"],
               xlabel="Время", ylabel="Популяция",
               title="Лотка-Вольтерра: $(parameter_names[idx])",
               linewidth=2,
               legend=:topright,
               size=(900, 400))

    hline!(plt, [fp_analysis.x_star], color=:green, linestyle=:dash, label="x* (теор.)")
    hline!(plt, [fp_analysis.y_star], color=:red, linestyle=:dash, label="y* (теор.)")

    savefig(plt, plotsdir(script_name, "lv_set$(idx)_dynamics.png"))
end

# ## СРАВНИТЕЛЬНЫЙ АНАЛИЗ ВСЕХ НАБОРОВ ПАРАМЕТРОВ

println("\n" * "="^80)
println(" " * "СРАВНИТЕЛЬНЫЙ АНАЛИЗ НАБОРОВ ПАРАМЕТРОВ")
println("="^80)

# Создаем сводную таблицу
comparison_df = DataFrame(
    Набор = Int[],
    Описание = String[],
    α = Float64[],
    β = Float64[],
    δ = Float64[],
    γ = Float64[],
    x_стационар = Float64[],
    y_стационар = Float64[],
    T_теор = Float64[],
    T_эксп = Float64[],
    x_мин = Float64[],
    x_макс = Float64[],
    x_среднее = Float64[],
    y_мин = Float64[],
    y_макс = Float64[],
    y_среднее = Float64[],
    отклонение_T = Float64[]
)

for (idx, r) in enumerate(results)
    push!(comparison_df, [
        idx,
        r.name,
        r.params[1],
        r.params[2],
        r.params[3],
        r.params[4],
        r.fp.x_star,
        r.fp.y_star,
        r.fp.period,
        1/r.dom_freq.prey,
        r.stats.min_prey,
        r.stats.max_prey,
        r.stats.mean_prey,
        r.stats.min_pred,
        r.stats.max_pred,
        r.stats.mean_pred,
        abs(1/r.dom_freq.prey - r.fp.period) / r.fp.period * 100
    ])
end

println("\n📋 Сводная таблица результатов:")
println(comparison_df)

# Сохранение сводной таблицы
CSV.write(datadir(script_name, "comparison_results.csv"), comparison_df)


# ## ВИЗУАЛИЗАЦИЯ СРАВНЕНИЯ ВСЕХ НАБОРОВ

# График 1: Сравнение периодов колебаний
plt_periods = bar(1:length(results),
                  [r.fp.period for r in results],
                  label="Теоретический период",
                  fillalpha=0.5,
                  xlabel="Номер набора параметров",
                  ylabel="Период колебаний",
                  title="Сравнение периодов колебаний")

bar!(plt_periods, 1:length(results),
     [1/r.dom_freq.prey for r in results],
     label="Экспериментальный период (жертвы)",
     fillalpha=0.5)

savefig(plt_periods, plotsdir(script_name, "comparison_periods.png"))

# График 2: Сравнение стационарных точек
x_stars = [r.fp.x_star for r in results]
y_stars = [r.fp.y_star for r in results]

plt_fixed = scatter(x_stars, y_stars,
                    label="Стационарные точки",
                    xlabel=L"x* (стационар жертв)",
                    ylabel=L"y* (стационар хищников)",
                    title="Стационарные точки для разных наборов параметров",
                    markersize=8,
                    markercolor=:blue,
                    grid=true)

for (idx, (x, y)) in enumerate(zip(x_stars, y_stars))
    annotate!(plt_fixed, x, y, text("  $idx", :left, 10))
end

savefig(plt_fixed, plotsdir(script_name, "comparison_fixed_points.png"))


println("\n✅ Моделирование успешно завершено!")
println("📁 Результаты сохранены в: $(datadir(script_name))")
println("📁 Графики сохранены в: $(plotsdir(script_name))")