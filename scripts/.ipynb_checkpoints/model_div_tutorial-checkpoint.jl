include("model_eqs_tutorial.jl")

gm_divC = @reaction_network gm begin 

    vimp(et, vt, s0, Kt),                                                                       ∅ => si
    nucat(em, vm, si, Km),                                                                      si => ∅
    ns * nucat(em, vm, si, Km),                                                                 ∅ => a

    ttrate(a, cr, ct, cm, cq, gmax, Kgamma),                                                    a => ∅

    # transcription 
    tx(we, thetax, a),                                                                      ∅ => mm
    tx(we, thetax, a),                                                                      ∅ => mt
    tx_q(wq, thetax, a, q=q, Kq=Kq, nq=nq),                                                     ∅ => mq
    tx(wr, thetax, a),                                                                      ∅ => mr

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
    tlr(a, nx, cm, gmax, Kgamma),                                                             cm => r + em + mm
    tlr(a, nx, ct, gmax, Kgamma),                                                             ct => r + et + mt
    tlr(a, nx, cq, gmax, Kgamma),                                                             cq => r + q + mq
    tlr(a, nr, cr, gmax, Kgamma),                                                             cr => r + mr + r

    # dilution
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * mt,                                              mt => ∅   
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * mm,                                             mm => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * mq,                                              mq => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * mr,                                             mr => ∅     
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * ct,                                              ct => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * cm,                                              cm => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * cq,                                              cq => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * cr,                                             cr => ∅    
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * et,                                              et => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * em,                                              em => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * q,                                               q => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * r,                                               r => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * si,                                              si => ∅
    lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) * a,                                               a => ∅

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


