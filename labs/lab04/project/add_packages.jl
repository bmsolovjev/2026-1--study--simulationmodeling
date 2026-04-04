#!/usr/bin/env julia
# add_packages.jl

using Pkg
Pkg.activate(".")

packages = [
    "DifferentialEquations",
    "Plots",
    "DataFrames",
    "CSV",
    "JLD2",
    "Literate",
    "IJulia",
    "BenchmarkTools",
    "Quarto",
    "Agents",
    "Random",
    "CairoMakie",
    "StatsBase",
    "Distributions",
    "BlackBoxOptim",
    "Statistics"
]

println("Установка базовых пакетов...")
Pkg.add(packages)

println("\n✅ Все пакеты установлены!")
