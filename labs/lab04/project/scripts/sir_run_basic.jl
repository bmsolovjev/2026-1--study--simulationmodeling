# # Модель SIR
# Загружаем проект и необходимые пакеты.
using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots
using JLD2

# Подключаем файл с определением модели.
include(srcdir("sir_model.jl"))

# ## Параметры эксперимента
# Задаём параметры для симуляции.
params = Dict(
    :Ns => [1000, 1000, 1000],
    :β_und => [0.5, 0.5, 0.5],
    :β_det => [0.05, 0.05, 0.05],
    :infection_period => 14,
    :detection_time => 7,
    :death_rate => 0.02,
    :reinfection_probability => 0.1,
    :Is => [0, 0, 1],
    :seed => 42,
    :n_steps => 100,
)

# ## Инициализация модели
# Создаём начальную модель с заданными параметрами.
model = initialize_sir(; params...)

# ## Запуск симуляции
# Подготавливаем массивы для хранения данных.
times = Int[]
S_vals = Int[]
I_vals = Int[]
R_vals = Int[]
total_vals = Int[]

# Запускаем симуляцию вручную на указанное количество шагов.
for step = 1:params[:n_steps]
    Agents.step!(model, 1)

    push!(times, step)
    push!(S_vals, susceptible_count(model))
    push!(I_vals, infected_count(model))
    push!(R_vals, recovered_count(model))
    push!(total_vals, total_count(model))
end

# ## Сохранение результатов
# Создаём DataFrame для удобного хранения.
agent_df = DataFrame(time = times, susceptible = S_vals, infected = I_vals, recovered = R_vals)
model_df = DataFrame(time = times, total = total_vals)

# Сохраняем данные в формате JLD2.
@save datadir("sir_basic_agent.jld2") agent_df
@save datadir("sir_basic_model.jld2") model_df

# ## Визуализация
# Строим график динамики эпидемии.
plot(agent_df.time, agent_df.susceptible,
     label = "Восприимчивые",
     xlabel = "Дни", ylabel = "Количество")
plot!(agent_df.time, agent_df.infected, label = "Инфицированные")
plot!(agent_df.time, agent_df.recovered, label = "Выздоровевшие")
plot!(agent_df.time, model_df.total, label = "Всего (включая умерших)", linestyle = :dash)

# Сохраняем график.
savefig(plotsdir("sir_basic_dynamics.png"))