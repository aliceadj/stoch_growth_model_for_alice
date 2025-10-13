function find_project_root(start_path=@__DIR__)
    current = abspath(start_path)
    while current != dirname(current)  # Stop at filesystem root
        if isfile(joinpath(current, "Project.toml"))
            return current
        end
        current = dirname(current)
    end
    error("Could not find project root (Project.toml not found)")
end
project_root = find_project_root()

import Pkg
Pkg.activate(joinpath(project_root, "stoch_growth_model_server"))
Pkg.instantiate()

using Distributed
nworkers1 = min(Sys.CPU_THREADS - 1, 31)  # Leave one core for system, max 31 workers - change this if there are other jobs needed to be run on the server 
addprocs(nworkers1)

@everywhere using Dates, Catalyst, JumpProcesses, Random, Distributions, DataFrames, DifferentialEquations, Arrow

@everywhere project_root = $project_root

@everywhere include(joinpath(project_root, "src/setup_funcs.jl"))
@everywhere include(joinpath(project_root, "src/models/catalyst_model/model_div.jl"))
@everywhere include(joinpath(project_root, "src/models/catalyst_model/stoch_callbacks.jl"))


abx_vec = [0,2,4,8,12]
ns_vec = [0.7481584182562523, 1.7623942844141665, 0.45586113648724746, 1.0365785284612912] 

param_combinations = [(ns, abx) for ns in ns_vec for abx in abx_vec]

date = Dates.format(Dates.now(), "ddmmyy_HHMM")

log_file = joinpath(project_root, "simulation_data/$date/simulation.log")
mkpath(dirname(log_file))

prob = DiscreteProblem(gm_divC, getX0(gm_divC, ss=true, parametrization="NatComms"), (0.0,1e5), getPars("molecs", "NatComms")[1:22])
jump_prob = JumpProblem(gm_divC, prob, Direct(), save_positions=(false, false), )

@everywhere jump_prob_template = $jump_prob

println("Starting simulations at $(date)")
@sync @distributed for i in eachindex(param_combinations)

    Random.seed!(i)

    state = SimState([],[],[])
    callbacks = make_callbacks(state)

    ns_val, abx_val = param_combinations[i]

    new_params = getPars("molecs", "NatComms"; ns=ns_val, abx=abx_val)[1:22]
    new_prob = remake(jump_prob_template, p=new_params)

    start_time = Dates.now()
    println("Starting run for ns=$(ns_val), abx=$(abx_val) at $(Dates.format(start_time, "HH:MM:SS"))")
    run_time = @elapsed sol = solve(new_prob, callback=CallbackSet(callbacks.fork_cb, callbacks.division_cb, callbacks.cellcycle_cb), saveat=0.1)
    end_time = Dates.now()

    df = DataFrame(sol)

    
    log_entry = "$(Dates.format(end_time, "yyyy-mm-ddTHH:MM:SS.sss")): abx=$(abx_val), ns=$(ns_val), time=$(round(run_time/60, digits=2)) minutes, $(round(run_time/60/60, digits=2)) hours"
    
    open(log_file, "a") do f
        println(f, log_entry)
    end

    savepath = joinpath(project_root, "simulation_data/$date/ns$(ns_val)")
    mkpath(savepath)  
    Arrow.write(joinpath(savepath, "seed$(i)_abx$(abx_val).arrow"), df, compress=:lz4)

    println("Finished run for ns=$(ns_val), abx=$(abx_val) at $(Dates.format(end_time, "HH:MM:SS"))")
end
println("Finished simulations at $(Dates.format(Dates.now(), "HH:MM:SS"))")


