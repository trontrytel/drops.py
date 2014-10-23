import inspect
from drops_py.output import output_lgr

class defaults:
  # common
  T  = None
  p  = None
  RH = None
  w  = None
  dt = None
  nt = None

  # lgrngn
  cloud_r_min, cloud_r_max = dict(zip(
    reversed(inspect.getargspec(output_lgr.__init__)[0]), 
    reversed(inspect.getargspec(output_lgr.__init__)[3]))
  )["cloud_rng"][0:2]
 
  sd_conc = None
  kappa   = None
  gstdv   = None
  meanr   = None
  n_tot   = None
  chem_SO2  = 0
  chem_O3   = 0
  chem_H2O2 = 0
