import libcloudphxx as libcl
import numpy as np
import h5py

class rhs_lgrngn:

  # ctor
  def __init__(self, dt, sd_conc, dry_distros, chem_gas = None, cloud_rng = None):
    opts_init = libcl.lgrngn.opts_init_t()
    opts_init.sd_conc_mean = sd_conc
    opts_init.dry_distros = dry_distros
    opts_init.dt = dt

    # TODO: not here
    self.cloud_rng = cloud_rng

    backend = libcl.lgrngn.backend_t.serial #TODO: as an option
    self.prtcls = libcl.lgrngn.factory(backend, opts_init)
    self.opts = libcl.lgrngn.opts_t()

    # turning off sedimentation and coalescence
    self.opts.sedi = False
    self.opts.coal = False

    # disabling substepping
    self.sstp_cond = 1

    if chem_gas != None:
      self.opts.chem = True
      self.opts.chem_gas = chem_gas

    # TODO: what's below should not be here...
    self.dt = dt


  # t=0 stuff
  def init(self, rhod, th_d, r_v):
    self.prtcls.init(th_d, r_v, rhod)

  def step(self, rhod, th_d, r_v, dot_th, dot_rv):
    th_d_copy = th_d.copy()
    r_v_copy = r_v.copy()
    self.prtcls.step_sync(self.opts, th_d_copy, r_v_copy, rhod)
    self.prtcls.step_async(self.opts)
    dot_th += (th_d_copy - th_d) / self.dt
    dot_rv += (r_v_copy - r_v) / self.dt

