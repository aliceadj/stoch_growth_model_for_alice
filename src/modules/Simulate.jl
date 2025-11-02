module Simulate

export define_jump_prob

include("StochModel.jl")
using .StochModel: growth_model

include("Setup.jl")
using .Setup: getPars, getX0

include("DetModel.jl")
using .DetModel: det_growth_model

using DifferentialEquations

function define_jump_prob(; units="molecs", parameterization="NatComms", ns=0.5, abx=0.0, ss=true, tspan=(0.0,1e5))

    gm_divC = growth_model()

    pars = getPars(units, parameterization; ns=ns, abx=abx)
    init = getX0(gm_divC, ss=ss, parametrization=parameterization)

    prob = DiscreteProblem(gm_divC, init, tspan, pars[1:22])
    jump_prob = JumpProblem(gm_divC, prob, Direct(), save_positions=(false, false), )

    return jump_prob
end

function simulate_ODE(; units="molecs", parameterization="NatComms", abx=0.0, ns=0.5)
    gm_dil = det_growth_model()
    pars = getPars(units, parameterization; abx=abx, ns=ns)
    init_dil = getX0(gm_dil)
    tspan = (0.1,1e9)

    prob_dil = ODEProblem(gm_dil, init_dil, tspan, pars)
    sol = solve(prob_dil, Rodas5())
    return sol 
end


end