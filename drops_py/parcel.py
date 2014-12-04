import libcloudphxx.common as libcom 
import analytical_equations as eq
import numpy as np

# p_d, th_d, r_v should contain initial values
#                and are overwritten!
def parcel(p_d, th_d, r_v, w, nt, outfreq, out, rhs, stats=None):

  # t=0 stuff
  rhod = eq.rhod_fun(p_d, th_d)
  rhs.init(rhod, th_d, r_v)

  # saving initial values
  #TODO: this call should be microphysics-independent - but now it relies on presence of prtcls in rhs!
  out.diag(rhs.prtcls, rhod, th_d, r_v, 0, stats=stats) 

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

    # recording values that were input for microphysics
    rhod_in = rhod.copy()
    th_d_in = th_d.copy()
    r_v_in  = r_v.copy()

    # applying the rhs
    th_d += rhs.dt * dot_th
    r_v  += rhs.dt * dot_rv
    rhod = eq.rhod_fun(p_d, th_d)

    # doing diagnostics / output
    out.diag(
      rhs.prtcls, 
      # values before adjustment by microphysics
      rhod_in, th_d_in, r_v_in,
      # values after adjustment by microphysics
      rhod, th_d, r_v, 
      (it+1) * rhs.dt, 
      stats = stats,
      save = (it+1) % outfreq == 0 # storing only every outfreq timesteps, 
                                   # calling anyway for housekeeping (e.g. S_max)
    ) 
