import libcloudphxx as libcl
import numpy as np
import h5py

class rhs_lgrngn:

  # ctor
  def __init__(self, dt, sd_conc, dry_distros, chem_gas = None):
    self.opts_init = libcl.lgrngn.opts_init_t()
    self.opts_init.sd_conc_mean = sd_conc
    self.opts_init.dry_distros = dry_distros
    self.opts_init.dt = dt
    self.dt = dt # as it is used from parcel.py :(

    self.opts = libcl.lgrngn.opts_t()

    # turning off sedimentation and coalescence
    self.opts.sedi = False
    self.opts.coal = False

    # disabling substepping
    self.opts_init.sstp_cond = 1

    if chem_gas != None:
      self.opts.chem = True
      self.opts.chem_gas = chem_gas


  # t=0 stuff
  def init(self, rhod, th_d, r_v):
    backend = libcl.lgrngn.backend_t.serial #TODO: as an option
    self.opts_init.RH_max = 1 - 1e-8 
    self.prtcls = libcl.lgrngn.factory(backend, self.opts_init)
    self.prtcls.init(th_d, r_v, rhod)

  def step(self, rhod, th_d, r_v, dot_th, dot_rv):
    th_d_copy = th_d.copy()
    r_v_copy = r_v.copy()
    self.prtcls.step_sync(self.opts, th_d_copy, r_v_copy, rhod)
    self.prtcls.step_async(self.opts)
    dot_th += (th_d_copy - th_d) / self.opts_init.dt
    dot_rv += (r_v_copy - r_v) / self.opts_init.dt

