include(joinpath(project_root, "src/models/model_eqs_tutorial.jl"))

gm_divC = @reaction_network gm begin 

    vimp(et, vt, s0, Kt),                                                                         ∅ => si
    vcat(em, vm, si, Km, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et),              si => ∅
    ns*vcat(em, si, vm, Km),           ∅ => a
    ttrate(a, cr, ct, cm, cq, gmax, Kgamma),        a => ∅ßß

    # transcription 
    tx_mtr(we, thetax, a),          ∅ => mm
    tx_mtr(we, thetax, a),          ∅ => mt
    tx_q(we, thetax, a, q=q, Kq=Kq, nq=nq),    ∅ => mq
    tx_mtr(we, thetax, a),          ∅ => mr

    # ribosome binding and unbinding
    ribo_bind(kb, mm, r),              mm + r => rmm
    ribo_bind(kb, mt, r),              mt + r => rmt
    ribo_bind(kb, mq, r),              mq + r => rmq
    ribo_bind(kb, mr, r),              mr + r => rmr

    ribo_unbind(ku, cm),                                                                         rmm => mm + r
    ribo_unbind(ku, ct),                                                                         rmt => mt + r
    ribo_unbind(ku, cq),                                                                         rmq => mq + r
    ribo_unbind(ku, cr),                                                                         rmr => mr + r

    # translation
    tlr_em(a, nx, cm, gmax, Kgamma),                                 rmm => r + em + mm
    tlr_et(a, nx, ct, gmax, Kgamma),        rmt => r + et + mt
    tlr_q(a, nx, cq, gmax, Kgamma),         rmq => r + q + mq
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


