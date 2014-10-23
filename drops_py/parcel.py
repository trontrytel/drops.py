import libcloudphxx.common as libcom 
import analytical_equations as eq
import numpy as np

# p_d, th_d, r_v should contain initial values
#                and are overwritten!
def parcel(p_d, th_d, r_v, w, nt, outfreq, out, rhs):

  # t=0 stuff
  rhod = eq.rhod_fun(p_d, th_d)
  rhs.init(rhod, th_d, r_v)

  # saving initial values
  out.diag(rhs.prtcls, rhod, th_d, r_v, 0)

  # Euler-like integration
  for it in range(nt):
    #TODO: update process name :)

    # first, adjusting thr pressure using hydrostatic law
    p_d += rhs.dt * (-libcom.g * rhod * w)

    # computing rhs for th and rv
    dot_th = np.array([0.])
    dot_rv = np.array([0.])
    rhod = eq.rhod_fun(p_d, th_d)
    rhs.step(rhod, th_d, r_v, dot_th, dot_rv)

    # applying the rhs
    th_d += rhs.dt * dot_th
    r_v  += rhs.dt * dot_rv
    rhod = eq.rhod_fun(p_d, th_d)

    # doing diagnostics / output
    if ((it+1)  % outfreq == 0):
      out.diag(rhs.prtcls, rhod, th_d, r_v, (it+1) * rhs.dt) 

