include(joinpath(project_root, "src/models/model_eqs.jl"))

gm_divC = @reaction_network gm begin 

    vimp(et, vt, s0, Kt),                                                                         ∅ => si
    nucat(em, vm, si, Km, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et),              si => ∅
    ns*nucat(em, vm, si, Km, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et),           ∅ => a
    ttrate(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),        a => ∅

    # transcription 
    tx_mtr(we, a, thetax, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),          ∅ => mm
    tx_mtr(we, a, thetax, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),          ∅ => mt
    tx_q(wq, a, thetax, Kq, nq, q, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, et, em),    ∅ => mq
    tx_mtr(wr, a, thetar, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),          ∅ => mr

    # ribosome binding and unbinding
    ribo_bind(kb, mm, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),              mm + r => rmm
    ribo_bind(kb, mt, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),              mt + r => rmt
    ribo_bind(kb, mq, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),              mq + r => rmq
    ribo_bind(kb, mr, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),              mr + r => rmr

    ribo_unbind(ku, rmm),                                                                         rmm => mm + r
    ribo_unbind(ku, rmt),                                                                         rmt => mt + r
    ribo_unbind(ku, rmq),                                                                         rmq => mq + r
    ribo_unbind(ku, rmr),                                                                         rmr => mr + r

    # translation
    tlr_em(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),        rmm => r + em + mm
    tlr_et(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),        rmt => r + et + mt
    tlr_q(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),         rmq => r + q + mq
    tlr_r(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em),         rmr => r + r + mr

    # degradation
    dm * mm,                                                                                      mm => ∅
    dm * mt,                                                                                      mt => ∅
    dm * mq,                                                                                      mq => ∅
    dm * mr,                                                                                      mr => ∅

    # antibiotic binding
    zm(rmm, abx, kon),                                                                            rmm => zmm
    zm(rmt, abx, kon),                                                                            rmt => zmt
    zm(rmq, abx, kon),                                                                            rmq => zmq
    zm(rmr, abx, kon),                                                                            rmr => zmr


end

# render(latexify(rn; form=:ode))

# gm_divC = complete(gm_divC)


# odesys_div = convert(ODESystem, gm_divC)
# odesys_div = complete(odesys_div)
# reactions(gm_divC)

# # division via callback 
# function condition(u, p, t) #condition for testing
#     # println("division!")
#     # calculateMassLive(u) - 7e7*2
#     mass(7549, u[7], u[8], u[9], u[10], u[11], 300, u[14], u[13], u[12]) - 7e7*2
# end

# function affect!(integrator) # affect testing
#     for i in 1:14
#         integrator.u[i] /= 2
#     end
#     # reset_aggregated_jumps!(integrator)
# end

# cb = ContinuousCallback(condition, affect!)

# division via jump 
# function mass_rate(u, p, t)
#     mass(7549.0, u[12], u[8], u[7], u[6], u[5], 300, u[9], u[10], u[11]) - 1e6 < 1e-4 ? 1.0 : 0.0 
# end
# function mass_affect!(integrator)
#     println("division!")
#     integrator.u.u ./= 2
# end

# division_jump = VariableRateJump(mass_rate, mass_affect!)
