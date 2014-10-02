from libcloudphxx.common import R_d, c_pd, g, p_1000
import numpy as np
import output

#TODO: should be moved to new file with physical equations etc. ??
# perfect gas for for dry air
def T_fun(p_d, th_d):
  return th_d * pow(p_d / p_1000, R_d / c_pd)                                                           
def rhod_fun(p_d, th_d):
  return p_d / R_d / T_fun(p_d, th_d)


# p_d, th_d, r_v should contain initial values
#                and are overwritten!
def parcel(p_d, th_d, r_v, w, nt, outfreq, rhs):

  # t=0 stuff
  rhod = rhod_fun(p_d, th_d)
  rhs.init(rhod, th_d, r_v)

  # preparing output
  it_out_ar = np.arange(0, nt+0.0001, outfreq) #TODO better (how?)
  out_file = output.output_lgr(rhs.outdir, it_out_ar)
  out_file.initial_info()
  # saving initial values
  out_file.diag(rhs.prtcls, rhod, th_d, r_v, 0, 0)

  # placing a quick-look gnuplot file in the output directory
  import os, shutil
  shutil.copyfile(os.path.dirname(__file__) + '/quicklook.gpi', rhs.outdir + '/quicklook.gpi')
  
  it_out = 1 # TODO better? 
  # Euler-like integration
  for it in range(nt):
    #TODO: update process name :)

    # first, adjusting thr pressure using hydrostatic law
    p_d += rhs.dt * (-g * rhod * w)

    # computing rhs for th and rv
    dot_th = np.array([0.])
    dot_rv = np.array([0.])
    rhod = rhod_fun(p_d, th_d)
    rhs.step(rhod, th_d, r_v, dot_th, dot_rv)

    # applying the rhs
    th_d += rhs.dt * dot_th
    r_v  += rhs.dt * dot_rv
    rhod = rhod_fun(p_d, th_d)

    # doing diagnostics / output
    if ((it+1) in it_out_ar):
      out_file.diag(rhs.prtcls, rhod, th_d, r_v, (it+1) * rhs.dt, it_out) 
      it_out += 1
