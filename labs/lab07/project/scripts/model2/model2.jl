using DrWatson
@quickactivate "project"

using ResumableFunctions
using ConcurrentSim
using Distributions
using Random
using StableRNGs
using DataFrames, CSV, Plots

RUNS = 5
N = 10
S = 3
NUM_REPAIRMEN = 1
SEED = 150
LAMBDA = 100
MU = 1

const rng = StableRNG(42)
const F = Exponential(LAMBDA)
const G = Exponential(MU)

global working_count = N + S
global busy_repairmen = 0
global waiting_for_repair = 0
working_machines = Tuple{Float64, Int}[]
repair_load = Tuple{Float64, Int}[]
queue_length = Tuple{Float64, Int}[]

@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    while true
        try
            @yield timeout(env, Inf)
        catch
        end
        @yield timeout(env, rand(rng, F))

        global working_count
        working_count -= 1
        push!(working_machines, (now(env), working_count))

        get_spare = take!(spares)
        @yield get_spare | timeout(env)
        if state(get_spare) != ConcurrentSim.idle
            @yield interrupt(value(get_spare))

            push!(working_machines, (now(env), working_count))
        else
            push!(working_machines, (now(env), working_count))
            throw(StopSimulation("No more spares!"))
        end

        global waiting_for_repair
        waiting_for_repair += 1
        push!(queue_length, (now(env), waiting_for_repair))

        @yield request(repair_facility)

        global waiting_for_repair
        waiting_for_repair -= 1

        global busy_repairmen
        busy_repairmen += 1
        push!(repair_load, (now(env), busy_repairmen))
        push!(queue_length, (now(env), waiting_for_repair))

        @yield timeout(env, rand(rng, G))
        @yield unlock(repair_facility)

        global busy_repairmen
        busy_repairmen -= 1
        push!(repair_load, (now(env), busy_repairmen))

        @yield put!(spares, active_process(env))

        global working_count
        working_count += 1
        push!(working_machines, (now(env), working_count))
    end
end

@resumable function start_sim(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    for i = 1:N
        proc = @process machine(env, repair_facility, spares)
        @yield interrupt(proc)
    end
    for i = 1:S
        proc = @process machine(env, repair_facility, spares)
        @yield put!(spares, proc)
    end
end

function sim_repair(; n::Int = N, s::Int = S, r::Int = NUM_REPAIRMEN)
    global working_machines = Tuple{Float64, Int}[]
    global repair_load = Tuple{Float64, Int}[]
    global queue_length = Tuple{Float64, Int}[]
    global working_count = n + s    # исправных = основные + запасные
    global busy_repairmen = 0
    global waiting_for_repair = 0

    push!(working_machines, (0.0, n + s))
    push!(repair_load, (0.0, 0))
    push!(queue_length, (0.0, 0))

    sim = Simulation()
    repair_facility = Resource(sim, r)
    spares = Store{Process}(sim)
    @process start_sim(sim, repair_facility, spares)
    msg = run(sim)
    stop_time = now(sim)
    return stop_time, copy(working_machines), copy(repair_load), copy(queue_length)
end

configs = [(N=10, S=3), (N=10, S=5), (N=15, S=3), (N=15, S=5)]
results = DataFrame(N = Int[], S = Int[], crash_time = Float64[])

for (n, s) in configs
    times = Float64[]
    for _ = 1:RUNS
        t_sim, _, _, _ = sim_repair(n = n, s = s)
        push!(times, t_sim)
    end
    push!(results, (n, s, mean(times)))
    println("N=$n, S=$s: среднее время до отказа = $(round(mean(times), digits=1)) часов")
end

CSV.write(datadir("scan_results.csv"), results)

t_sim, working_log, load_log, queue_log = sim_repair()

plot_limit = min(100.0, t_sim)

times_w = [x[1] for x in working_log]
machines = [x[2] for x in working_log]

idx_w = findall(t -> t <= plot_limit, times_w)
times_w_zoom = times_w[idx_w]
machines_zoom = machines[idx_w]

if length(times_w_zoom) > 0
    p1 = plot(times_w_zoom, machines_zoom,
              title = "Число исправных машин (первые $plot_limit часов)",
              xlabel = "Время (часы)", ylabel = "Исправные машины",
              label = "Исправные машины", lw = 2, marker = :circle, markersize = 3,
              ylims = (N + S - 2, N + S + 1))
    hline!([N + S], linestyle = :dash, color = :blue, label = "N+S = $(N+S) (всего)")
    hline!([N], linestyle = :dash, color = :green, label = "N = $N (основные)")
    savefig(plotsdir("working_machines.png"))
end

times_l = [x[1] for x in load_log]
loads = [x[2] for x in load_log]

idx_l = findall(t -> t <= plot_limit, times_l)
times_l_zoom = times_l[idx_l]
loads_zoom = loads[idx_l]

if length(times_l_zoom) > 0
    p2 = plot(times_l_zoom, loads_zoom,
              title = "Загрузка ремонтника (первые $plot_limit часов)",
              xlabel = "Время (часы)", ylabel = "Занятые ремонтники",
              label = "Занят ремонтник", lw = 2, marker = :circle, markersize = 3,
              ylims = (0, NUM_REPAIRMEN + 0.5))
    hline!([NUM_REPAIRMEN], linestyle = :dash, color = :red, label = "Всего ремонтников")
    savefig(plotsdir("repair_load.png"))
end

times_q = [x[1] for x in queue_log]
queues = [x[2] for x in queue_log]

idx_q = findall(t -> t <= plot_limit, times_q)
times_q_zoom = times_q[idx_q]
queues_zoom = queues[idx_q]

if length(times_q_zoom) > 0 && maximum(queues_zoom) > 0
    p3 = plot(times_q_zoom, queues_zoom,
              title = "Длина очереди на ремонт (первые $plot_limit часов)",
              xlabel = "Время (часы)", ylabel = "Машин в очереди",
              label = "Очередь", lw = 2, marker = :circle, markersize = 3,
              ylims = (0, max(maximum(queues_zoom) + 0.5, 1)))
    savefig(plotsdir("queue_length.png"))
else
    p3 = plot([0], [0],
              title = "Длина очереди на ремонт (очередь пуста)",
              xlabel = "Время (часы)", ylabel = "Машин в очереди",
              label = "", lw = 2, marker = :circle,
              ylims = (0, 1))
    savefig(plotsdir("queue_length.png"))
end

println("\n--- Статистика за первые $plot_limit часов ---")
if length(times_w_zoom) > 0
    println("Число событий изменения числа исправных машин: $(length(times_w_zoom))")
    println("Минимум исправных машин: $(minimum(machines_zoom))")
    println("Максимум исправных машин: $(maximum(machines_zoom))")
end
if length(times_l_zoom) > 0
    println("Событий загрузки ремонтника: $(length(times_l_zoom))")
end
if length(times_q_zoom) > 0
    println("Событий в очереди: $(length(times_q_zoom))")
end

all_loads = [x[2] for x in load_log]
all_queues = [x[2] for x in queue_log]
avg_load = mean(all_loads)
avg_queue = mean(all_queues)
println("\n--- Мониторинг (весь прогон) ---")
println("Общее время симуляции: $(round(t_sim, digits=1)) часов")
println("Средняя загрузка ремонтников: $(round(avg_load, digits=2)) из $NUM_REPAIRMEN")
println("Средняя длина очереди на ремонт: $(round(avg_queue, digits=2)) машин")
println("Коэффициент загрузки: $(round(avg_load / NUM_REPAIRMEN * 100, digits=1))%")

lambda_sys = N / LAMBDA
mu_repair = 1 / MU

println("\n" * "="^60)
println("АНАЛИТИЧЕСКОЕ СРАВНЕНИЕ")
println("="^60)
println("N (основных машин) = $N")
println("S (запасных машин) = $S")
println("Ремонтников = $NUM_REPAIRMEN")
println("λ одной машины = $(1/LAMBDA) отказов/час")
println("λ системы = $lambda_sys отказов/час")
println("μ ремонтника = $mu_repair ремонтов/час")
println("ρ = λ/μ = $(round(lambda_sys/mu_repair, digits=3))")

λ = lambda_sys
μ = mu_repair
r = λ / μ

sum_pi = sum([r^k for k = 0:S])
pi_S = r^S / sum_pi
failure_rate = λ * pi_S
time_to_failure = 1 / failure_rate

println("\nВероятности состояний:")
for k = 0:S
    pi = r^k / sum_pi
    println("  $k неисправных машин: $(round(pi * 100, digits=2))%")
end

println("\nИнтенсивность отказа системы: $(round(failure_rate, digits=6)) отказов/час")
println("Аналитическое среднее время до отказа: $(round(time_to_failure, digits=1)) часов")

if r < 1
    Lq = r^2 / (1 - r)
    println("Аналитическая средняя очередь (M/M/1): $(round(Lq, digits=3)) машин")
    println("Аналитическая загрузка ремонтника: $(round(r * 100, digits=1))%")
end

sim_times = Float64[]
for _ = 1:RUNS
    local t, _, _, _ = sim_repair(n = N, s = S, r = NUM_REPAIRMEN)
    push!(sim_times, t)
end
sim_avg = mean(sim_times)
println("\n--- Сравнение ---")
println("Симуляция (среднее по $RUNS прогонам): $(round(sim_avg, digits=1)) часов")
println("Относительная погрешность: $(round(abs(sim_avg - time_to_failure) / time_to_failure * 100, digits=1))%")

println("\nРезультаты сохранены в data/ и plots/")
