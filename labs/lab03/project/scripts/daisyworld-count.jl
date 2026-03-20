# ## Анализ динамики популяции маргариток
# Загружаем проект и необходимые пакеты
using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, CairoMakie

# ## Подключение модели
include(srcdir("daisyworld.jl"))

# ## Определение функций для сбора данных
# Функции для подсчёта чёрных и белых маргариток
black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]

# ## Запуск модели с изменённой солнечной активностью
# Модель запускается на 1000 шагов при солнечной светимости 1.0
model = daisyworld(; solar_luminosity=1.0)
agent_df, model_df = run!(model, 1000; adata)

# ## Визуализация динамики популяций
# Создаём график изменения численности маргариток во времени
figure = Figure(size=(600, 400))
ax = figure[1, 1] = Axis(figure, xlabel="tick", ylabel="daisy count")

blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black], color=:black)
whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white], color=:orange)

Legend(figure[1, 2], [blackl, whitel], ["black", "white"], labelsize=12)
save(plotsdir("daisy_count.png"), figure)