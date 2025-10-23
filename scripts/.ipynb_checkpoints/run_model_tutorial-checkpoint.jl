using Catalyst, JumpProcesses, DifferentialEquations
using DataFrames, CSV, Plots, Dates


function find_project_root(start_path=@__DIR__)
    current = abspath(start_path)
    while current != dirname(current)
        if isfile(joinpath(current, "Project.toml"))
            return current
        end
        current = dirname(current)
    end
    error("Could not find project root (Project.toml not found)")
end

project_root = find_project_root()

include(joinpath(project_root, "src/models/model_eqs_tutorial.jl"))
include(joinpath(project_root, "src/models/catalyst_model/model_div_tutorial.jl"))
include(joinpath(project_root, "src/setup_funcs.jl"))

params = getPars("molecs", "NatComms"; abx = 0.0)
u0 = Dict(getX0(gm_divC; ss = false, parametrization = "NatComms"))

for s in (:ct, :cm, :cq, :cr)
    u0[s] = get(u0, s, 0.0)
end

println("Starting stochastic simulation")

dprob = DiscreteProblem(gm_divC, u0, (0.0, 1_000.0), params)
jprob = JumpProblem(gm_divC, dprob, Direct(); save_positions = (false, false))

sol = solve(jprob, SSAStepper(); saveat = 1.0)

println("Simulation complete")

df = DataFrame(sol)
CSV.write(joinpath(project_root, "run_model_tutorial.csv"), df)

plot(sol;
     vars = [:a, :si, :r, :em, :et, :q],
     lw = 1.5,
     xlabel = "Time",
     ylabel = "Molecule count",
     legend = :outerright,
     title = "Stochastic Simulation of gm_divC")

savefig(joinpath(project_root, "run_model_tutorial.png"))
println("Saved outputs")