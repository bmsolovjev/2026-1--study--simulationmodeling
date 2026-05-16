using DrWatson
@quickactivate "project"

using StableRNGs
using Distributions
using ConcurrentSim
using ResumableFunctions
using DataFrames, CSV, Plots


# Параметры модели

rng = StableRNG(123)
num_customers = 10       # общее число клиентов
num_servers = 2          # количество серверов
mu = 1.0 / 2             # интенсивность обслуживания
lam = 0.9                # интенсивность прибытия
arrival_dist = Exponential(1 / lam)    # распределение интервалов прибытия
service_dist = Exponential(1 / mu)     # распределение времени обслуживания


# Хранилище для сбора статистики

log = DataFrame(
    customer_id = Int[],
    event = String[],
    time = Float64[]
)

# Словари для хранения времён каждого клиента
arrival_times = Dict{Int, Float64}()
enter_times = Dict{Int, Float64}()
exit_times = Dict{Int, Float64}()


# Поведение клиента

@resumable function customer(
    env::Environment,
    server::Resource,
    id::Integer,
    t_a::Float64,
    d_s::Distribution,
)
    @yield timeout(env, t_a)                        # прибытие клиента
    println("Customer $id arrived: ", now(env))
    push!(log, (id, "arrival", now(env)))
    arrival_times[id] = now(env)
    
    @yield request(server)                          # начало обслуживания
    println("Customer $id entered service: ", now(env))
    push!(log, (id, "enter_service", now(env)))
    enter_times[id] = now(env)
    
    service_time = rand(rng, d_s)
    @yield timeout(env, service_time)               # обслуживание
    
    @yield unlock(server)                           # выход из системы
    println("Customer $id exited service: ", now(env))
    push!(log, (id, "exit_service", now(env)))
    exit_times[id] = now(env)
end

# Запуск симуляции

function setup_and_run()
    global log = DataFrame(customer_id = Int[], event = String[], time = Float64[])
    global arrival_times = Dict{Int, Float64}()
    global enter_times = Dict{Int, Float64}()
    global exit_times = Dict{Int, Float64}()
    
    sim = Simulation()
    server = Resource(sim, num_servers)
    arrival_time = 0.0
    for i = 1:num_customers
        arrival_time += rand(rng, arrival_dist)
        @process customer(sim, server, i, arrival_time, service_dist)
    end
    run(sim)
    return log
end


# Выполнение

data = setup_and_run()
println("Моделирование завершено. Собрано ", nrow(data), " событий.")


CSV.write(datadir("queue_log.csv"), data)


# График 1: Временная шкала обслуживания клиентов

p1 = plot(legend = :topleft, title = "Временная шкала обслуживания", 
          xlabel = "Время", ylabel = "Клиент")
for id in 1:num_customers
    if haskey(arrival_times, id) && haskey(enter_times, id) && haskey(exit_times, id)
        arr = arrival_times[id]
        ent = enter_times[id]
        ext = exit_times[id]
        
        plot!(p1, [arr, ent], [id, id], lw = 2, color = :orange, label = nothing)
        plot!(p1, [ent, ext], [id, id], lw = 4, color = :steelblue, label = nothing)
        scatter!(p1, [arr], [id], marker = :circle, color = :orange, label = nothing)
        scatter!(p1, [ext], [id], marker = :circle, color = :red, label = nothing)
    end
end
scatter!(p1, [NaN], [NaN], marker = :circle, color = :orange, label = "Ожидание")
scatter!(p1, [NaN], [NaN], marker = :circle, color = :steelblue, label = "Обслуживание")
scatter!(p1, [NaN], [NaN], marker = :circle, color = :red, label = "Выход")
savefig(plotsdir("queue_timeline.png"))


# График 2: Длительность ожидания и обслуживания

wait_times = Float64[]
service_times = Float64[]
for id in 1:num_customers
    if haskey(arrival_times, id) && haskey(enter_times, id) && haskey(exit_times, id)
        push!(wait_times, enter_times[id] - arrival_times[id])
        push!(service_times, exit_times[id] - enter_times[id])
    end
end

p2 = bar(
    [wait_times service_times],
    label = ["Ожидание" "Обслуживание"],
    title = "Время ожидания и обслуживания по клиентам",
    xlabel = "Клиент",
    ylabel = "Время",
    legend = :topleft,
    color = [:orange :steelblue]
)

savefig(plotsdir("queue_wait_service.png"))
