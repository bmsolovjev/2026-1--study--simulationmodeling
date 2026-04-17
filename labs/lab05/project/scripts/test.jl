# # Моделирование задачи об обедающих философах
#
# Демонстрация стохастического поведения двух сетей Петри:
# 1. Классическая сеть (подвержена deadlock)
# 2. Сеть с арбитром (свободна от deadlock)

# ## Подготовка окружения

using DrWatson
@quickactivate "project"
include(srcdir("DiningPhilosophers.jl"))
using .DiningPhilosophers
using DataFrames, CSV, Plots

# ## Параметры моделирования

N = 5           # Количество философов
tmax = 50.0     # Время моделирования

# ---
# ## 1. Классическая сеть (без арбитра)
#
# В классической постановке каждый философ захватывает левую вилку, 
# затем правую. При совпадении моментов захвата возможен deadlock.

println("=== Классическая сеть (без арбитра) ===")
net_classic, u0_classic, _ = build_classical_network(N)
df_classic = simulate_stochastic(net_classic, u0_classic, tmax)

# Сохранение результатов
CSV.write(datadir("dining_classic.csv"), df_classic)

# Проверка на deadlock
dead = detect_deadlock(df_classic, net_classic)
println("Deadlock обнаружен: $dead")

# Визуализация эволюции маркировки
plot_classic = plot_marking_evolution(df_classic, N)
savefig(plotsdir("classic_simulation.png"))

# ---
# ## 2. Сеть с арбитром
#
# Арбитр управляет очерёдностью захвата вилок, исключая 
# возникновение циклического ожидания.

println("\n=== Сеть с арбитром ===")
net_arb, u0_arb, _ = build_arbiter_network(N)
df_arb = simulate_stochastic(net_arb, u0_arb, tmax)

# Сохранение результатов
CSV.write(datadir("dining_arbiter.csv"), df_arb)

# Проверка на deadlock (ожидается false)
dead_arb = detect_deadlock(df_arb, net_arb)
println("Deadlock обнаружен: $dead_arb")

# Визуализация эволюции маркировки
plot_arb = plot_marking_evolution(df_arb, N)
savefig(plotsdir("arbiter_simulation.png"))

# ---
# ## Выводы
#
# - Классическая сеть демонстрирует возможность deadlock при 
#   стохастическом моделировании.
# - Сеть с арбитром гарантированно избегает deadlock.
# - Результаты сохранены в `data/` и `plots/`.
#
# Для повторного анализа загрузите данные командой:
# ```julia
# df = CSV.read(datadir("dining_classic.csv"), DataFrame)
# ```