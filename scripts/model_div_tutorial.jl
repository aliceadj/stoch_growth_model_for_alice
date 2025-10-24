include("model_eqs_tutorial.jl")

gm_divC = @reaction_network gm begin 

    vimp(et, vt, s0, Kt),                                                                       ∅ => si
    nucat(em, vm, si, Km), si => ∅
    ns * nucat(em, vm, si, Km), ∅ => a

    ttrate(a, cr, ct, cm, cq, gmax, Kgamma),                                                    a => ∅

    # transcription 
    tx_mtr(we, thetax, a),                                                                      ∅ => mm
    tx_mtr(we, thetax, a),                                                                      ∅ => mt
    tx_q(we, thetax, a, q=q, Kq=Kq, nq=nq),                                                     ∅ => mq
    tx_mtr(we, thetax, a),                                                                      ∅ => mr

    # ribosome binding and unbinding
    ribo_bind(kb, mm, r),                                                                       mm + r => cm
    ribo_bind(kb, mt, r),                                                                       mt + r => ct
    ribo_bind(kb, mq, r),                                                                       mq + r => cq
    ribo_bind(kb, mr, r),                                                                       mr + r => cr

    ribo_unbind(ku, cm),                                                                        cm => mm + r
    ribo_unbind(ku, ct),                                                                        ct => mt + r
    ribo_unbind(ku, cq),                                                                        cq => mq + r
    ribo_unbind(ku, cr),                                                                        cr => mr + r

    # translation
    tlr_em(a, nx, cm, gmax, Kgamma),                                                            cm => r + em + mm
    tlr_et(a, nx, ct, gmax, Kgamma),                                                            ct => r + et + mt
    tlr_q(a, nx, cq, gmax, Kgamma),                                                             cq => r + q + mq
    tlr_r(gmax, a, Kgamma, nr, r, cm, ct, cq, cr, zmm, zmt, zmq, zmr, nx, q, et, em),       cr => r + r + mr

    # dilution
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              mt => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                             mm => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              mq => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                             mr => ∅     
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              ct => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              cm => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              cq => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                             cr => ∅    
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              et => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              em => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                               q => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                               r => ∅

    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                              si => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M),                                               a => ∅

    # degradation
    dm * mm,                                                                                      mm => ∅
    dm * mt,                                                                                      mt => ∅
    dm * mq,                                                                                      mq => ∅
    dm * mr,                                                                                      mr => ∅

    # antibiotic binding
    zm(cm, abx, kon),                                                                            cm => zmm
    zm(ct, abx, kon),                                                                            ct => zmt
    zm(cq, abx, kon),                                                                            cq => zmq
    zm(cr, abx, kon),                                                                            cr => zmr
end


