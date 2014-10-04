from libcloudphxx.common import g
import analytical_equations as eq
import numpy as np
import output


# p_d, th_d, r_v should contain initial values
#                and are overwritten!
def parcel(p_d, th_d, r_v, w, nt, outfreq, rhs):

  # t=0 stuff
  rhod = eq.rhod_fun(p_d, th_d)
  rhs.init(rhod, th_d, r_v)

  # preparing output
  time_out = rhs.dt * np.arange(0, nt+1, outfreq) # nt+1 to include nt in the time_out
  out_file = output.output_lgr(rhs.outdir, time_out)
  # saving initial values
  out_file.diag(rhs.prtcls, rhod, th_d, r_v, 0)

  # placing a quick-look gnuplot file in the output directory
  import os, shutil
  shutil.copyfile(os.path.dirname(__file__) + '/quicklook.gpi', rhs.outdir + '/quicklook.gpi')
  
  # Euler-like integration
  for it in range(nt):
    #TODO: update process name :)

    # first, adjusting thr pressure using hydrostatic law
    p_d += rhs.dt * (-g * rhod * w)

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
      out_file.diag(rhs.prtcls, rhod, th_d, r_v, (it+1) * rhs.dt) 

