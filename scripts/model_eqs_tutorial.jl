v1(nr, r, cm, ct, cq, cr, zmm, zmt, zmq, zmr, nx, q, et, em)  = 1e8/mass(nr, r, cm, ct, cq, cr, zmm, zmt, zmq, zmr, nx, q, et, em)
mass(nr, r, cm, ct, cq, cr, zmm, zmt, zmq, zmr, nx, q, et, em)  = nr * (r + cm + ct + cq + cr + zmm + zmt + zmq + zmr) + nx * (q + et + em)

vimp(et, vt, s0, Kt)           = (et * vt * (s0) / ((Kt) + (s0))) 
nucat(em, vm, si, Km)           = (em * vm * (si) / (Km + si))
vcat_a(em, vm, si, Km, ns)     = (ns*(em * vm * si / ((Km) + si)))

tx_q(wq, theta, a, q, Kq, nq) = (wq * a / (theta + a)) / (1 + (q / Kq)^nq)
tx(w, theta, a) = w * a / (theta + a)


ribo_bind(kb, mx, r)   = (kb) * mx * r 
ribo_unbind(ku, cx)      = ku * cx 

# translation
gamma_fun(gmax, a, Kgamma) = gmax * a / (Kgamma + a)

# tlr_r(gmax, a, Kgamma, nr, r, cm, ct, cq, cr, zmm, zmt, zmq, zmr, nx, q, et, em) = (gamma_fun(gmax, a, Kgamma) / nx) * cr
tlr(a, nx, cx, gmax, Kgamma) = (gamma_fun(gmax, a, Kgamma) / nx) * cx

ttrate(a, cr, ct, cm, cq, gmax, Kgamma) = (cr + ct + cm + cq) * gamma_fun(gmax, a, Kgamma)

lam(a, cr, ct, cm, cq, gmax, Kgamma, M, nr, nx, r, zmm, zmt, zmq, zmr, q, et, em) = (ttrate(a, cr, ct, cm, cq, gmax, Kgamma) / 1e8) 



zm(c_x, abx, kon)                                                     = c_x * abx * kon  