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

ns_vec = [0.7481584182562523, 1.7623942844141665, 0.45586113648724746, 1.0365785284612912] # according to Christoph's data (see ns_estimation in src/models/catalyst_model)

# add workers so each simulation runs on a different process
using Distributed
nworkers1 = length(ns_vec)
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
    sol = simulate_ODE(; units="molecs", parameterization="NatComms", abx=0.0, ns=ns)
    init_conds[ns] = sol[end]
end

# load in jump problem
jump_prob = define_jump_prob(units="molecs", parameterization="NatComms", tspan=(0.0,1e5)) #1e5 or 1

@everywhere jump_prob_template = $jump_prob

println("Starting simulations at $(date)")
@sync @distributed for i in eachindex(ns_vec)
    # sets a seed to make simulations reproducible
    Random.seed!(i)

    state = SimState([],[],[])
    callbacks = make_callbacks(state)

    # get parameter values for this run
    ns_val = ns_vec[i]
    abx_val = 0.0

    new_params = getPars("molecs", "NatComms"; ns=ns_val, abx=abx_val)[1:22]
    # remake the jump problem with new parameters and initial conditions
    new_prob = remake(jump_prob_template, p=new_params, u0=init_conds[ns_val])

    start_time = Dates.now()
    println("Starting run for ns=$(ns_val), abx=0.0 at $(Dates.format(start_time, "HH:MM:SS"))")
    run_time = @elapsed sol = solve(new_prob, callback=CallbackSet(callbacks.fork_cb, callbacks.division_cb, callbacks.cellcycle_cb), saveat=10/60)
    end_time = Dates.now()

    df = DataFrame(sol)

    # record time it took for each simulation in a file 
    log_entry = "$(Dates.format(end_time, "yyyy-mm-ddTHH:MM:SS.sss")): abx=$(abx_val), ns=$(ns_val), time=$(round(run_time/60, digits=2)) minutes, $(round(run_time/60/60, digits=2)) hours"
    
    open(log_file, "a") do f
        println(f, log_entry)
    end

    # save the results in a folder
    savepath = joinpath(project_root, "simulation_data/$date")
    mkpath(savepath)  
    Arrow.write(joinpath(savepath, "seed$(i)_ns$(ns_val)_abx$(abx_val).arrow"), df, compress=:lz4)

    println("Finished run for ns=$(ns_val), abx=$(abx_val) at $(Dates.format(end_time, "HH:MM:SS"))")
end
println("Finished simulations at $(Dates.format(Dates.now(), "HH:MM:SS"))")


