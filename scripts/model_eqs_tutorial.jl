#v1(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)                          = 1e8/mass(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)
#mass(nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em)                        = nr * (r + rmm + rmt + rmq + rmr + zmm + zmt + zmq + zmr) + nx * (q + et + em)

vimp(et, vt, s0, Kt)           = (et * vt * (s0) / ((Kt) + (s0))) 
nucat(em, vm, si, Km)           = (em * vm * (si) / (Km + si))

#tx_mtrtx(w, theta, a; q=nothing, Kq=nothing, nq=nothing) = (w * a / (theta + a)) / (1 + (q / Kq)^nq)
tx_mtrtx(w, theta, a; q=0.0, Kq=1.0, nq=1.0) = (w * a / (theta + a)) / (1 + (q / Kq)^nq)
tx_mtr(w, theta, a) = tx_mtrtx(w, theta, a)
#tx_q(w, theta, a, q, Kq, nq) = (w * a / (theta + a)) / (1 + (q / (Kq ))^nq)
tx_q(w, theta, a; q=0.0, Kq=1.0, nq=1.0) = (w * a / (theta + a)) / (1 + (q / Kq)^nq)

ribo_bind(kb, mx, r)           = (kb) * mx * r 

ribo_unbind(ku, cx)      = ku * cx 

# translation
gamma_fun(gmax, a, Kgamma) = gmax * a / (Kgamma + a)
tlr(a, nx, cx, gmax, Kgamma) = (gamma_fun(gmax, a, Kgamma) / nx) * cx
tlr_em(a, nx, cx, gmax, Kgamma) = tlr(a, nx, cx, gmax, Kgamma)
tlr_et(a, nx, cx, gmax, Kgamma) = tlr(a, nx, cx, gmax, Kgamma)
tlr_q(a, nx, cx, gmax, Kgamma)  = tlr(a, nx, cx, gmax, Kgamma)
tlr_r(gmax, a, Kgamma, nr, r, rmm, rmt, rmq, rmr, zmm, zmt, zmq, zmr, nx, q, et, em) = (gamma_fun(gmax, a, Kgamma) / nx) * rmr


ttrate(a, cr, ct, cm, cq, gmax, Kgamma) = (cr + ct + cm + cq) * gamma_fun(gmax, a, Kgamma)
lam(a, cr, ct, cm, cq, gmax, Kgamma, M) = ttrate(a, cr, ct, cm, cq, gmax, Kgamma) / M


zm(c_x, abx, kon)                                                     = c_x * abx * kon  