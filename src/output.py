import os
import h5py
import numpy as np

class output:

  def __init__(self, outdir):
    self.outdir_hdf = outdir + "/hdf_output/"
    try:
      os.mkdir(outdir)
    except OSError:
      pass
    try:
      os.mkdir(self.outdir_hdf)
    except OSError:
      pass
    self.out_snd = open(outdir + "/sounding.txt", mode='w')
    self.out_dry = open(outdir + "/spec_dry.txt", mode='w')
    self.out_wet = open(outdir + "/spec_wet.txt", mode='w')
    self.bins_dry = 1e-6 * pow(10, -3 + np.arange(40) * .1)
    self.bins_wet = 1e-6 * pow(10, -3 + np.arange(25) * .2)

  def initial_info(self):
    self.out_snd.write(u"#rhod [kg/m3]\tth_d [K] (theta dry!)\tr_v [kg/kg] (mixing ratio)\tM0 [TODO]\tM1 [TODO]\tM2 [TODO]\tM3 [TODO]\tS_VI [kg/kg]\tH [kg/kg]\tSO2 [kg/kg]\n")
    self.out_dry.write(u"#r_d [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
    self.out_wet.write(u"#r_w [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")


  def diag(self, prtcls, rhod, th_d, r_v, time):
    #opening hdf files (can't be in init: new file for each diag(time)) TODO: should be changed?
    fhdf_spec = h5py.File(self.outdir_hdf + "spec_drywet_" + str(int(time)) + ".hdf", mode='w')
    fhdf_sound = h5py.File(self.outdir_hdf + "sounding_" + str(int(time)) + ".hdf", mode='w')
    # helper for gnuplot-readable output                                
    def save(out, xx, yy):
      for x, y in zip(xx, yy):
        out.write(u"%g\t%g\n" % (x, y))
      out.write(u"\n\n") # gnuplot treats "\n\n" as dataset separator (plot ... index n)
    # outputting spec_dry.txt
    bins = np.empty(self.bins_dry.size - 1)
    for i in range(0, bins.size) :
      prtcls.diag_dry_rng(self.bins_dry[i], self.bins_dry[i+1])
      prtcls.diag_dry_mom(0)
      bins[i] = np.frombuffer(prtcls.outbuf())
      save(self.out_dry, self.bins_dry[0:-1], bins)
    # outputting to hdf 
    dry_r_h5 = fhdf_spec.create_dataset("dry_radius", (bins.size,), dtype='f')
    dry_r_h5.attrs["Units"] = "m"
    dry_r_h5[...] = self.bins_dry[:-1]
    dry_n_h5 = fhdf_spec.create_dataset("dry_number", (bins.size,), dtype='f')
    dry_n_h5[...] = bins
    dry_n_h5.attrs["Units"] = "1/kg"

    # outputting spec_wet.txt 
    bins = np.empty(self.bins_wet.size - 1)
    for i in range(0, bins.size) :
      prtcls.diag_wet_rng(self.bins_wet[i], self.bins_wet[i+1])
      prtcls.diag_wet_mom(0)
      bins[i] = np.frombuffer(prtcls.outbuf())
    save(self.out_wet, self.bins_wet[0:-1], bins)
    # outputting to hdf
    wet_r_h5 = fhdf_spec.create_dataset("wet_radius", (bins.size,), dtype='f')
    wet_r_h5[...] = self.bins_wet[:-1]
    wet_r_h5.attrs["Units"] = "m"
    wet_n_h5 = fhdf_spec.create_dataset("wet_number", (bins.size,), dtype='f')
    wet_n_h5[...] = bins
    wet_n_h5.attrs["Units"] = "1/kg"



    # outputting sounding.txt
    self.out_snd.write(u"%g" % (rhod))
    self.out_snd.write(u"\t%g" % (th_d))
    self.out_snd.write(u"\t%g" % (r_v))
    # outputting to hdf 
    rho_h5 = fhdf_sound.create_dataset("rhod", (rhod.size,), dtype='f')
    rho_h5[...] = rhod
    rho_h5.attrs["Units"] = "kg/m3"
    thd_h5 = fhdf_sound.create_dataset("thd", (th_d.size,), dtype='f')
    thd_h5[...] = th_d
    thd_h5.attrs["Units"] = "K"
    rv_h5 = fhdf_sound.create_dataset("rv", (r_v.size,), dtype='f')
    rv_h5[...] = r_v
    rv_h5.attrs["Units"] = "kg/kg"


    ## cloud water
    prtcls.diag_wet_rng(.5e-6, 25e-6)
    for k in range(0,4):
      prtcls.diag_wet_mom(k)
      self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))


    ## chem stuff
    #prtcls.diag_wet_rng(0,1) # 0 ... 1 m #TODO: consider a select-all option?
    #prtcls.diag_chem(lgrngn.chem_species_t.S_VI)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
    #prtcls.diag_chem(lgrngn.chem_species_t.H)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
    #prtcls.diag_chem(lgrngn.chem_species_t.SO2)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))

    self.out_snd.write(u"\n")

