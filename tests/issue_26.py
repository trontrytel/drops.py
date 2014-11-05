#!/usr/bin/python

import tempfile, os, sys

outdir = tempfile.mkdtemp()
srcdir = sys.argv[1]

os.system(srcdir + "/drops.py --outdir " + outdir + 
  " --T 288 --p 100000 --RH 0.999 --w 0.62593 --outfreq 10000" + 
  " --dt 0.5 --nt 320 lgrngn --sd_conc 128 --kappa 0.6052" + 
  " --n_tot 930759410 317786280 --meanr 1.1429e-07 4.3822e-08" + 
  " --gstdv 1.4745 1.4466"
)

f = open(outdir + "/stats.gpi")
for line in f:
  flds = line.split(' = ')
  if flds[0] == 'S_max_RH':
    print flds[1]
    if float(flds[1]) < 1: 
      print "maximal supersaturation less than 1!"
      sys.exit(1)
    else: sys.exit(0)
print "no maximal supersaturation has been found!"
sys.exit(1)
