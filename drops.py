#!/usr/bin/python

from argparse import ArgumentParser
from drops_py import rhs_lgrngn, parcel
import libcloudphxx.common as libcom 
import libcloudphxx as libcl
import numpy as np
import math 


# TODO: move to a separate file
class defaults:
  chem_SO2  = 0
  chem_O3   = 0
  chem_H2O2 = 0
  T = None
  p = None
  # TODO...

# TODO: move to a separate file
# Ghan, S., Guzman, G., and Abdul-Razzak, H. 1998 Competition between sea salt
#   and sulfate particles as cloud condensation nuclei,
# J. Atmos. Sci., 55, 3340-3347
# doi:10.1175/1520-0469(1998)055<3340:CBSSAS>2.0.CO;2
class defaults_Ghan_et_al_1998(defaults):
  T = 280
  p = 100000
  # TODO...

# just a few constants not repeat them below
desc = 'drops.py - a parcel model based on libcloudph++'
chcs = ['Ghan_et_al_1998']
dhlp = 'default parameter set'

# defining command-line options parser - first without help
# so that it would not be empty when calling parse_known_args() below
prsr = ArgumentParser(add_help=False, description=desc)
prsr.add_argument('--defaults', choices=chcs, help=dhlp)

# first parsing the "defaults" only - to know the defaults
args = prsr.parse_known_args() # TODO: take care of prefix mathing rules = "def" will also match here :(
defclass = 'defaults'
if args[0].defaults is not None:
  defclass = "defaults_" + args[0].defaults
ctor = globals()[defclass]
dflts = ctor()

# redefining prsr, this time with help
prsr = ArgumentParser(add_help=True, description=desc)
prsr.add_argument('--defaults', choices=chcs, help=dhlp)

# one subparser per one microphysics scheme
sprsr = prsr.add_subparsers()

## common options
prsr.add_argument('--outdir', required=True, help='output directory')

## lgrngn options
prsr_lgr = sprsr.add_parser('lgrngn')
prsr_lgr.add_argument('--sd_conc',   type=float, required=True, help='number of super droplets')
prsr_lgr.add_argument('--T',         type=float, required=(dflts.T is None), default=dflts.T, help='initial temperature [K]')
prsr_lgr.add_argument('--p',         type=float, required=(dflts.p is None), default=dflts.p, help='initial pressure [Pa]')
prsr_lgr.add_argument('--RH',        type=float, required=True, help='initial relative humidity [1]')
prsr_lgr.add_argument('--w',         type=float, required=True, help='vertical velocity [m/s]')
prsr_lgr.add_argument('--dt',        type=float, required=True, help='timestep [s]')
prsr_lgr.add_argument('--nt',        type=int,   required=True, help='number of timesteps')
prsr_lgr.add_argument('--outfreq',   type=int,   required=True, help='output frequency (every outfreq timesteps)')
prsr_lgr.add_argument('--kappa',     type=float, required=True, help='aerosol hygroscopicity parameter [1]')
prsr_lgr.add_argument('--n_tot',     type=float, required=True, help='aerosol concentration @STP [m-3]')
prsr_lgr.add_argument('--meanr',     type=float, required=True, help='aerosol mean dry radius [m]')
prsr_lgr.add_argument('--gstdv',     type=float, required=True, help='aerosol geometric standard deviation [1]')
prsr_lgr.add_argument('--chem_SO2',  type=float, default=dflts.chem_SO2,  help='SO2 volume concentration [1]')
prsr_lgr.add_argument('--chem_O3',   type=float, default=dflts.chem_O3,   help='O3 volume concentration [1]')
prsr_lgr.add_argument('--chem_H2O2', type=float, default=dflts.chem_H2O2, help='H2O2 volume concentration [1]')

## blk_2m options
prsr_b2m = sprsr.add_parser('blk_2m')
#TODO...

args = prsr.parse_args()

# computing state variables
p_v = np.array([args.RH * libcom.p_vs(args.T)])
p_d = args.p - p_v
r_v = libcom.eps * p_v / p_d
th_d = args.T * pow(libcom.p_1000 / p_d, libcom.R_d / libcom.c_pd)

class lognormal:
  def __init__(self, n_tot, meanr, gstdv):
    self.meanr = meanr
    self.stdev = gstdv
    self.n_tot = n_tot
 
  def __call__(self, lnr):
    return self.n_tot * math.exp(
      -pow((lnr - math.log(self.meanr)), 2) / 2 / pow(math.log(self.stdev),2)
    ) / math.log(self.stdev) / math.sqrt(2*math.pi);

# performing the simulation
rhs = rhs_lgrngn.rhs_lgrngn(
  args.outdir, 
  args.dt, 
  args.sd_conc, 
  { 
    args.kappa : lognormal(args.n_tot, args.meanr, args.gstdv)
  }
#,
#  {
#    libcl.lgrngn.chem_species_t.SO2  : args.chem_SO2,
#    libcl.lgrngn.chem_species_t.O3   : args.chem_O3,
#    libcl.lgrngn.chem_species_t.H2O2 : args.chem_H2O2
#  }
)
parcel.parcel(p_d, th_d, r_v, args.w, args.nt, args.outfreq, rhs)

# outputting a setup.gpi file
out = open(args.outdir + '/setup.gpi', mode='w')
for key, val in vars(args).iteritems():
  if key != "outdir":
    out.write(u"%s = %g\n" % (key, float(val)))
