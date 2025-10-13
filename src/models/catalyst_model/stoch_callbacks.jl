mutable struct SimState
    fork_creation_times::Vector{Float64}
    fork_creation_times_save::Vector{Float64}
    cellcycle_count::Vector{Float64}
end


function make_callbacks(state::SimState)

    fork_rate = @inline function (u, t, integrator)
        n_ori = 2^length(state.fork_creation_times)
        current_mass = mass(7549, u[7], u[8], u[9], u[10], u[11], u[15], u[16], u[17], u[18], 300, u[14], u[13], u[12])
        return n_ori < (1e-9 * current_mass)
    end

    fork_affect! = @inline function (integrator)
        push!(state.fork_creation_times, integrator.t)
        push!(state.fork_creation_times_save, integrator.t)
        println("Fork created at time $(integrator.t)")
    end

    division_rate = @inline function (u, t, integrator)
        return any(abs(t - (t_fork + 60)) < 1e-3 for t_fork in state.fork_creation_times)
    end

    division_affect! = @inline function (integrator)
        pDiv = rand(Beta(200.0, 200.0))
        for i in 1:18
            n = integrator.u[i]
            integrator.u[i] = n > 0 ? float(rand(Binomial(round(Int, n), pDiv))) : 0
        end
        # Keep only forks not ready for division
        state.fork_creation_times = filter(t_fork -> integrator.t < t_fork + 60 - 1e-3, state.fork_creation_times)
        println("Cell divided at $(integrator.t)")
        push!(state.cellcycle_count, integrator.t)
        reset_aggregated_jumps!(integrator, update_jump_params=false)
    end

    cellcycle_rate = @inline (u, t, integrator) -> length(state.cellcycle_count) >= 500

    cellcycle_affect! = @inline function (integrator)
        terminate!(integrator)
        println("Terminated at $(integrator.t) after $(length(state.cellcycle_count)) divisions")
    end

    return (
        fork_cb = DiscreteCallback(fork_rate, fork_affect!),
        division_cb = DiscreteCallback(division_rate, division_affect!),
        cellcycle_cb = DiscreteCallback(cellcycle_rate, cellcycle_affect!)
    )
end

