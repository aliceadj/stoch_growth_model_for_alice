# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
using Catalyst
using DifferentialEquations
using DataFrames
using CSV
using Plots

project_root = dirname(@__DIR__)

include(joinpath(project_root, "src/models/model_eqs_tutorial.jl"))
include(joinpath(project_root, "src/models/catalyst_model/model_div_tutorial.jl"))
include(joinpath(project_root, "src/setup_funcs.jl"))  
params = getPars("uM", "NatComms"; abx=0.0)
u0 = getX0(gm_divC; ss=false, parametrization="NatComms")


println("Solving ODE system")
prob = ODEProblem(gm_divC, u0, (0.0, 1e3), params)
sol = solve(prob, Tsit5(); saveat=1.0)
println("Simulation complete")


df = DataFrame(sol)
CSV.write(joinpath(project_root, "model_div_results.csv"), df)

plt = plot(sol;
    vars=[:a, :si, :r, :em, :et, :q],
    lw=2,
    xlabel="Time",
    ylabel="Concentration (µM)",
    legend=:outerright,
    title="Model DIV ODE Simulation"
)
savefig(joinpath(project_root, "model_div_ODE_results.png"))
println("Saved outputs: model_div_results.csv and model_div_ODE_results.png")


# %%


# %%
