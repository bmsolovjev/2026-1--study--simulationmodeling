# ## Инициализация проекта
# Активируем проект DrWatson и загружаем необходимые пакеты
using DrWatson
@quickactivate "project"
using Agents, DataFrames, Plots, CairoMakie

# ## Загрузка модели
# Подключаем файл с определением модели daisyworld
include(srcdir("daisyworld.jl"))

# ## Создание и настройка модели
# Инициализируем модель и определяем цвет для отображения агентов
model = daisyworld()
daisycolor(a::Daisy) = a.breed

# ## Параметры визуализации
# Настраиваем внешний графиков
plotkwargs = (
    agent_color=daisycolor, agent_size=20, agent_marker='✿',
    heatarray=:temperature, heatkwargs=(colorrange=(-20, 60),)
)

# ## Визуализация на разных шагах
# Шаг 0: начальное состояние
plt1, _ = abmplot(model; plotkwargs...)
save(plotsdir("daisy_step001.png"), plt1)

# Шаг 5: первые изменения
step!(model, 5)
plt2, _ = abmplot(model; heatarray=model.temperature, plotkwargs...)
save(plotsdir("daisy_step005.png"), plt2)

# Шаг 45: итоговое состояние
step!(model, 40)
plt3, _ = abmplot(model; heatarray=model.temperature, plotkwargs...)
save(plotsdir("daisy_step045.png"), plt3)