include("src/modules/Setup.jl")
using .Setup: find_project_root

# get the root of the project so files are saved in the correct location
project_root = find_project_root()

# activate the appropriate Julia environment depending if working on own device (assumed mac) or the server (linux)
import Pkg
if Sys.KERNEL == :Linux
    Pkg.activate(joinpath(project_root, "stoch_growth_model_server"))
elseif Sys.KERNEL == :Darwin
    Pkg.activate(project_root)
end

abx_list = [2,4,8,12]
ns_vec = [0.7481584182562523, 1.7623942844141665, 0.45586113648724746, 1.0365785284612912] # according to Christoph's data (see ns_estimation in src/models/catalyst_model)

# add workers so each simulation runs on a different process
using Distributed
nworkers1 = length(ns_vec) * length(abx_list)
addprocs(nworkers1)

# load in things across all workers
@everywhere using Dates, Catalyst, JumpProcesses, Random, Distributions, DataFrames, DifferentialEquations, Arrow

@everywhere project_root = $project_root

# load in modules and functions that are needed 
@everywhere include(joinpath(project_root, "src/modules/Callbacks.jl"))
@everywhere using .Callbacks: SimState, make_callbacks

@everywhere include(joinpath(project_root, "src/modules/Setup.jl"))
@everywhere using .Setup: getPars

include(joinpath(project_root, "src/modules/Simulate.jl"))
using .Simulate: define_jump_prob, simulate_ODE

date = Dates.format(Dates.now(), "ddmmyy_HHMM")

log_file = joinpath(project_root, "simulation_data/$date/simulation.log")
mkpath(dirname(log_file))


# set initial conditions to be the steady states of the ODE simulation in the same conditions
init_conds = Dict{Float64,Vector{Float64}}()
for ns in ns_vec
    for abx in abx_list
        sol = simulate_ODE(; units="molecs",
                            parameterization="NatComms",
                            abx=abx,
                            ns=ns)
        init_conds_local[(ns, abx)] = sol[end]
    end
end

# load in jump problem
base_ns  = ns_vec[1]
base_abx = abx_list[1]

base_params = getPars("molecs", "NatComms"; ns=base_ns, abx=base_abx)[1:22]

jump_prob = define_jump_prob(units="molecs", parameterization="NatComms", tspan=(0.0,1)) #1e5
param_combinations = [(ns, abx) for ns in ns_vec for abx in abx_list]

@everywhere jump_prob_template = $jump_prob

println("Starting simulations at $(date)")
@sync @distributed for idx in eachindex(param_combinations)
    ns_val, abx_val = param_combinations[idx]

    # Reproducible randomness
    Random.seed!(idx)

    state = SimState([], [], [])
    callbacks = make_callbacks(state)

    # Parameters for this run
    new_params = getPars("molecs", "NatComms"; ns=ns_val, abx=abx_val)[1:22]
    u0 = init_conds[(ns_val, abx_val)]

    new_prob = remake(jump_prob_template, p=new_params, u0=u0)

    start_time = Dates.now()
    println("Starting run for ns=$(ns_val), abx=$(abx_val) at $(Dates.format(start_time, "HH:MM:SS"))")

    run_time = @elapsed begin
        sol = solve(new_prob,
                    callback=CallbackSet(callbacks.fork_cb,
                                         callbacks.division_cb,
                                         callbacks.cellcycle_cb),
                    saveat=10/60)
        # Save solution to Arrow
        df = DataFrame(sol)

        savepath = $save_root
        mkpath(savepath)

        Arrow.write(joinpath(savepath,
                             "seed$(idx)_ns$(ns_val)_abx$(abx_val).arrow"),
                    df, compress=:lz4)
    end

    end_time = Dates.now()

    # Log entry (note: multiple workers append; that's ok for simple logs)
    log_entry = "$(Dates.format(end_time, "yyyy-mm-ddTHH:MM:SS.sss")): " *
                "abx=$(abx_val), ns=$(ns_val), " *
                "time=$(round(run_time/60, digits=2)) minutes, " *
                "$(round(run_time/60/60, digits=2)) hours"

    open($log_file, "a") do f
        println(f, log_entry)
    end

    println("Finished run for ns=$(ns_val), abx=$(abx_val) at $(Dates.format(end_time, "HH:MM:SS"))")
end
println("Finished simulations at $(Dates.format(Dates.now(), "HH:MM:SS"))")


