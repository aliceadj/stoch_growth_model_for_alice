module DetModel

export det_growth_model

using Catalyst  

vimp(et, vt, s0, Kt)           = (et * vt * (s0) / ((Kt) + (s0))) 
nucat(em, vm, si, Km)          = (em * vm * (si) / ((Km) + (si))) 
vcat_a(em, vm, si, Km, ns)     = (ns*(em * vm * si / ((Km) + si))) 

tx_emt(we, a, thetax)          = (we ) * a / ((thetax ) + a)
tx_q(wq, a, thetax, Kq, nq, q) = ((wq ) * a / ((thetax ) + a)) / (1 + (q / (Kq )) ^ nq)
tx_r(wr, a, thetar)            = (wr ) * a / ((thetar ) + a)

ribo_bind_em(kb, mm, r)        = (kb) * mm * r  
ribo_bind_et(kb, mt, r)        = (kb) * mt * r  
ribo_bind_q(kb, mq, r)         = (kb) * mq * r  
ribo_bind_r(kb, mr, r)         = (kb) * mr * r  

ribo_unbind_em(ku, rmm)              = ku * rmm #*v(mass)
ribo_unbind_et(ku, rmt)              = ku * rmt #*v(mass)
ribo_unbind_q(ku, rmq)               = ku * rmq  #*v(mass)
ribo_unbind_r(ku, rmr)               = ku * rmr  #*v(mass)

# translation
gamma(gmax, a, Kgamma)            = (gmax) * a / ((Kgamma ) + a)
tlr_em(gmax, a, Kgamma, nx, rmm)  = ((gamma(gmax, a, Kgamma) / nx) * rmm) #*v(mass)
tlr_et(gmax, a, Kgamma, nx, rmt)  = ((gamma(gmax, a, Kgamma) / nx) * rmt) #*v(mass)
tlr_q(gmax, a, Kgamma, nx, rmq)   = ((gamma(gmax, a, Kgamma) / nx) * rmq)  #*v(mass)
tlr_r(gmax, a, Kgamma, nr, rmr)   = ((gamma(gmax, a, Kgamma) / nr) * rmr)  #*v(mass)


ttrate(gmax, a, Kgamma, rmr, rmt, rmm, rmq) = ((rmr + rmt + rmm + rmq)*gamma(gmax, a, Kgamma))
lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, mass)    = ttrate(gmax, a, Kgamma, rmr, rmt, rmm, rmq)/mass

zm(c_x, abx, kon) = c_x * abx * kon  

function det_growth_model()
    gm_dilC = @reaction_network gm begin 

        vimp(et, vt, s0, Kt),                                    ∅ => si
        nucat(em, vm, si, Km),                                   si => ∅
        vcat_a(em, vm, si, Km, ns),                              ∅ => a
        ttrate(gmax, a, Kgamma, rmr, rmt, rmm, rmq),             a => ∅

        # transcription 
        tx_emt(we, a, thetax),                                   ∅ => mm
        tx_emt(we, a, thetax),                                   ∅ => mt
        tx_q(wq, a, thetax, Kq, nq, q),                          ∅ => mq
        tx_r(wr, a, thetar),                                     ∅ => mr

        # ribosome binding and unbinding
        ribo_bind_em(kb, mm, r),                                 mm + r => rmm
        ribo_bind_et(kb, mt, r),                                 mt + r => rmt
        ribo_bind_q(kb, mq, r),                                  mq + r => rmq
        ribo_bind_r(kb, mr, r),                                  mr + r => rmr

        ribo_unbind_em(ku, rmm),                                 rmm => mm + r
        ribo_unbind_et(ku, rmt),                                 rmt => mt + r
        ribo_unbind_q(ku, rmq),                                  rmq => mq + r
        ribo_unbind_r(ku, rmr),                                  rmr => mr + r

        # translation
        tlr_em(gmax, a, Kgamma, nx, rmm),                        rmm => r + em + mm
        tlr_et(gmax, a, Kgamma, nx, rmt),                        rmt => r + et + mt
        tlr_q(gmax, a, Kgamma, nx, rmq),                         rmq => r + q + mq
        tlr_r(gmax, a, Kgamma, nr, rmr),                         rmr => r + r + mr

        # degradation
        dm * mm,                                                  mm => ∅
        dm * mt,                                                  mt => ∅
        dm * mq,                                                  mq => ∅
        dm * mr,                                                  mr => ∅

        # antibiotic binding
        zm(rmm, abx, kon),                                        rmm => zmm
        zm(rmt, abx, kon),                                        rmt => zmt
        zm(rmq, abx, kon),                                        rmq => zmq
        zm(rmr, abx, kon),                                        rmr => zmr
        
        
        # dilution
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*a,         a => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*si,        si => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*mm,        mm => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*mt,        mt => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*mq,        mq => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*mr,        mr => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*rmm,       rmm => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*rmt,       rmt => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*rmq,       rmq => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*rmr,       rmr => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*em,        em => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*et,        et => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*q,         q => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*r,         r => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*zmm,       zmm => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*zmt,       zmt => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*zmq,       zmq => ∅
        lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, Mref)*zmr,       zmr => ∅


    end
    return gm_dilC
end
end

