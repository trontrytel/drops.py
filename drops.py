#!/usr/bin/python

from argparse import ArgumentParser
from drops_py import parcel, rhs_lgrngn
from libcloudphxx.common import p_vs, eps, p_1000, R_d, c_pd
from libcloudphxx.lgrngn import chem_species_t
from numpy import array as arr_t
from math import exp, log, sqrt, pi

# command-line options
prsr = ArgumentParser(description='drops.py - a parcel model based on libcloudph++')
sprsr = prsr.add_subparsers()

## common options
prsr.add_argument('--outdir', required=True, help='output directory')

## lgrngn options
prsr_lgr = sprsr.add_parser('lgrngn')
prsr_lgr.add_argument('--sd_conc',   type=float, required=True, help='number of super droplets')
prsr_lgr.add_argument('--T',         type=float, required=True, help='initial temperature [K]')
prsr_lgr.add_argument('--p',         type=float, required=True, help='initial pressure [Pa]')
prsr_lgr.add_argument('--RH',        type=float, required=True, help='initial relative humidity [1]')
prsr_lgr.add_argument('--w',         type=float, required=True, help='vertical velocity [m/s]')
prsr_lgr.add_argument('--dt',        type=float, required=True, help='timestep [s]')
prsr_lgr.add_argument('--nt',        type=int,   required=True, help='number of timesteps')
prsr_lgr.add_argument('--outfreq',   type=int,   required=True, help='output frequency (every outfreq timesteps)')
prsr_lgr.add_argument('--kappa',     type=float, required=True, help='aerosol hygroscopicity parameter [1]')
prsr_lgr.add_argument('--n_tot',     type=float, required=True, help='aerosol concentration @STP [m-3]')
prsr_lgr.add_argument('--meanr',     type=float, required=True, help='aerosol mean dry radius [m]')
prsr_lgr.add_argument('--gstdv',     type=float, required=True, help='aerosol geometric standard deviation [1]')
prsr_lgr.add_argument('--chem_SO2',  type=float, default=0,     help='SO2 volume concentration [1]')
prsr_lgr.add_argument('--chem_O3',   type=float, default=0,     help='O3 volume concentration [1]')
prsr_lgr.add_argument('--chem_H2O2', type=float, default=0,     help='H2O2 volume concentration [1]')

## blk_2m options
prsr_b2m = sprsr.add_parser('blk_2m')
#TODO...

args = prsr.parse_args()

# computing state variables
p_v = arr_t([args.RH * p_vs(args.T)])
p_d = args.p - p_v
r_v = eps * p_v / p_d
th_d = args.T * pow(p_1000 / p_d, R_d / c_pd)

class lognormal:
  def __init__(self, n_tot, meanr, gstdv):
    self.meanr = meanr
    self.stdev = gstdv
    self.n_tot = n_tot
 
  def __call__(self, lnr):
    return self.n_tot * exp(
      -pow((lnr - log(self.meanr)), 2) / 2 / pow(log(self.stdev),2)
    ) / log(self.stdev) / sqrt(2*pi);

# outputting a setup.gpi file
out = open(args.outdir + '/setup.gpi', mode='w')
for key, val in vars(args).iteritems():
  if key != "outdir":
    out.write(u"%s = %g\n" % (key, float(val)))

# performing the simulation
rhs = rhs_lgrngn.rhs_lgrngn(
  args.outdir, 
  args.dt, 
  args.sd_conc, 
  { 
    args.kappa : lognormal(args.n_tot, args.meanr, args.gstdv)
  },
  {
    chem_species_t.SO2  : args.chem_SO2,
    chem_species_t.O3   : args.chem_O3,
    chem_species_t.H2O2 : args.chem_H2O2
  }
)
parcel.parcel(p_d, th_d, r_v, args.w, args.nt, args.outfreq, rhs)
