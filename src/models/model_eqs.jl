
# v1(mass) = 1e8/mass
# mass(nr, r, rmm, rmt, rmq, rmr, nx, q, et, em) = nr * (r + rmm + rmt + rmq + rmr) + nx * (q + et + em)

# vimp(et, vt, s0, Kt, mass)           = (et * vt * (s0) / ((Kt) + (s0))) 
# nucat(em, vm, si, Km, mass)          = (em * vm * (si) / ((Km)/v1(mass) + (si))) 
# vcat_a(em, vm, si, Km, ns, mass)     = (ns*(em * vm * si / ((Km) + si))) 

# tx_emt(we, a, thetax, mass)          = (we /v1(mass)) * a / ((thetax /v1(mass)) + a)
# tx_q(wq, a, thetax, Kq, nq, q, mass) = ((wq /v1(mass)) * a / ((thetax /v1(mass)) + a)) / (1 + (q / (Kq /v1(mass))) ^ nq)
# tx_r(wr, a, thetar, mass)            = (wr /v1(mass)) * a / ((thetar /v1(mass)) + a)

# ribo_bind_em(kb, mm, r, mass)        = (kb) * mm * r * v1(mass) 
# ribo_bind_et(kb, mt, r, mass)        = (kb) * mt * r * v1(mass) 
# ribo_bind_q(kb, mq, r, mass)         = (kb) * mq * r * v1(mass) 
# ribo_bind_r(kb, mr, r, mass)         = (kb) * mr * r * v1(mass) 

# ribo_unbind_em(ku, rmm)              = ku * rmm #*v(mass)
# ribo_unbind_et(ku, rmt)              = ku * rmt #*v(mass)
# ribo_unbind_q(ku, rmq)               = ku * rmq  #*v(mass)
# ribo_unbind_r(ku, rmr)               = ku * rmr  #*v(mass)

# # translation
# gamma(gmax, a, Kgamma, mass)            = (gmax) * a / ((Kgamma /v1(mass)) + a)
# tlr_em(gmax, a, Kgamma, nx, rmm, mass)  = ((gamma(gmax, a, Kgamma, mass) / nx) * rmm) #*v(mass)
# tlr_et(gmax, a, Kgamma, nx, rmt, mass)  = ((gamma(gmax, a, Kgamma, mass) / nx) * rmt) #*v(mass)
# tlr_q(gmax, a, Kgamma, nx, rmq, mass)   = ((gamma(gmax, a, Kgamma, mass) / nx) * rmq)  #*v(mass)
# tlr_r(gmax, a, Kgamma, nr, rmr, mass)   = ((gamma(gmax, a, Kgamma, mass) / nr) * rmr)  #*v(mass)


# ttrate(gmax, a, Kgamma, rmr, rmt, rmm, rmq, mass) = ((rmr + rmt + rmm + rmq)*gamma(gmax, a, Kgamma, mass))
# lam(gmax, a, Kgamma, rmr, rmt, rmm, rmq, mass)    = ttrate(gmax, a, Kgamma, rmr, rmt, rmm, rmq, mass)/mass





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

# zm_diss(z_x, koff)                                                    = koff * z_x