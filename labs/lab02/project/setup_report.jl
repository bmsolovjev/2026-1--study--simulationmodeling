# Универсальная настройка для отчётов
using DrWatson

# Определяем путь к проекту (поднимаемся из report в labXX/project)
function find_project_path()
    current_dir = pwd()
    if basename(current_dir) == "report"
        return joinpath(dirname(current_dir), "project")
    else
        # fallback для других случаев
        return joinpath(dirname(dirname(current_dir)), "project")
    end
end

project_path = find_project_path()
println("Активируем проект: $project_path")
quickactivate(project_path)

# Устанавливаем имя скрипта для отчёта
script_name = "report"
mkpath(plotsdir(script_name))
mkpath(datadir(script_name))

println("Проект активирован: ", projectname())