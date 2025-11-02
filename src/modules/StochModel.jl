module StochModel

export growth_model, mass

using Catalyst 

v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)                          = 1e8/mass(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)
mass(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)                        = nr * (r + rmm + rmt + rmq + rmr + zmm + zmt + zmq + zmr) + nx * (q + et + em)

vimp(et, vt, s0, Kt)           = (et * vt * (s0) / ((Kt) + (s0))) 
nucat(em, vm, si, Km, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et)           = (em * vm * (si) / ((Km)/v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) + (si))) 

tx_mtr(wx, a, thetax, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)       = (wx /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)) * a / ((thetax /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)) + a)
tx_q(wq, a, thetax, Kq, nq, q, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, et, em) = ((wq /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)) * a / ((thetax /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)) + a)) / (1 + (q / (Kq /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em))) ^ nq)

ribo_bind(kb, mx, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)           = (kb) * mx * r * v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) 

ribo_unbind(ku, rmx)                                                  = ku * rmx 

# translation
gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)      = (gmax) * a / ((Kgamma /v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)) + a)
tlr_em(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)     = ((gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) / nx) * rmm) 
tlr_et(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)     = ((gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) / nx) * rmt) 
tlr_q(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)      = ((gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) / nx) * rmq) 
tlr_r(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)      = ((gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) / nr) * rmr) 


ttrate(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)     = ((rmr + rmt + rmm + rmq)*gamma(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em))
lam(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)        = ttrate(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)/mass(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)

zm(c_x, abx, kon)                                                     = c_x * abx * kon  

function growth_model()


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

    return gm_divC
end

end