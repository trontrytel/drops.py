import os
import h5py
import numpy as np

class output_lgr:

  def __init__(self, outdir, time_out):
    self.outdir_hdf = outdir + "/hdf_output/"
    try:
      os.mkdir(outdir)
    except OSError:
      pass
    try:
      os.mkdir(self.outdir_hdf)
    except OSError:
      pass
    self.time = time_out
    self.hdf_spec = h5py.File(self.outdir_hdf + "spec_drywet.hdf", mode='w')
    self.hdf_sound = h5py.File(self.outdir_hdf + "sounding.hdf", mode='w')
    self.out_snd = open(outdir + "/sounding.txt", mode='w')
    self.out_dry = open(outdir + "/spec_dry.txt", mode='w')
    self.out_wet = open(outdir + "/spec_wet.txt", mode='w')
    self.bins_dry = 1e-6 * pow(10, -3 + np.arange(40) * .1)
    self.bins_wet = 1e-6 * pow(10, -3 + np.arange(25) * .2)
    # creating hdf varibales  
    # concentratio per size of dry particles
    # TODO: this concentration can change in time???
    self.dry_n_h5 = self.hdf_spec.create_dataset("dry_number", 
                                 (self.time.size, self.bins_dry.size - 1), dtype='f')
    # concentratio per size of wet particles   
    self.wet_n_h5 = self.hdf_spec.create_dataset("wet_number", 
                                 (self.time.size, self.bins_wet.size - 1), dtype='f') 
    # sounding: dry density, dry pot. temp., water vapour mix. rat.
    self.rho_h5 = self.hdf_sound.create_dataset("rhod", (self.time.size,), dtype='f') 
    self.thd_h5 = self.hdf_sound.create_dataset("thd", (self.time.size,), dtype='f')
    self.rv_h5 = self.hdf_sound.create_dataset("rv", (self.time.size,), dtype='f')    


  def initial_info(self): # TODO should be part of __init__? 
    self.out_snd.write(u"#rhod [kg/m3]\tth_d [K] (theta dry!)\tr_v [kg/kg] (mixing ratio)\tM0 [TODO]\tM1 [TODO]\tM2 [TODO]\tM3 [TODO]\tS_VI [kg/kg]\tH [kg/kg]\tSO2 [kg/kg]\n")
    self.out_dry.write(u"#r_d [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
    self.out_wet.write(u"#r_w [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")

    # attaching scales, units etc. to the hdf variables 
    self.hdf_spec["time"] = self.time
    self.hdf_spec["bins_dry"] = self.bins_dry[:-1] # TODO should be 0.5 * (bins_wet[:-1] + bins_wet[1:]) ??    
    self.hdf_spec["bins_wet"] = self.bins_wet[:-1]

    self.dry_n_h5.dims.create_scale(self.hdf_spec["time"], "time")
    self.dry_n_h5.dims[0].attach_scale(self.hdf_spec["time"])
    self.dry_n_h5.dims[0].label = 's'
    self.dry_n_h5.dims.create_scale(self.hdf_spec["bins_dry"], "dry_radius")
    self.dry_n_h5.dims[1].attach_scale(self.hdf_spec["bins_dry"])
    self.dry_n_h5.dims[1].label = 'm'
    self.dry_n_h5.attrs["Units"] = "1/kg"

    self.wet_n_h5.dims.create_scale(self.hdf_spec["time"], "time")
    self.wet_n_h5.dims[0].attach_scale(self.hdf_spec["time"])
    self.wet_n_h5.dims[0].label = 's'
    self.wet_n_h5.dims.create_scale(self.hdf_spec["bins_wet"], "wet_radius")
    self.wet_n_h5.dims[0].attach_scale(self.hdf_spec["bins_wet"])
    self.wet_n_h5.dims[0].label = 'm'
    self.wet_n_h5.attrs["Units"] = "1/kg"

    self.rho_h5.dims.create_scale(self.hdf_spec["time"], "time")
    self.rho_h5.dims[0].attach_scale(self.hdf_spec["time"])
    self.rho_h5.dims[0].label = 's'
    self.rho_h5.attrs["Units"] = "kg/m3"

    self.thd_h5.dims.create_scale(self.hdf_spec["time"], "time")
    self.thd_h5.dims[0].attach_scale(self.hdf_spec["time"])
    self.thd_h5.dims[0].label = 's'
    self.thd_h5.attrs["Units"] = "K"

    self.rv_h5.dims.create_scale(self.hdf_spec["time"], "time")
    self.rv_h5.dims[0].attach_scale(self.hdf_spec["time"])
    self.rv_h5.dims[0].label = 's'
    self.rv_h5.attrs["Units"] = "kg/kg"


  def diag(self, prtcls, rhod, th_d, r_v, time, it_out):

    def save(out, xx, yy):
      for x, y in zip(xx, yy):
        out.write(u"%g\t%g\n" % (x, y))
      out.write(u"\n\n") # gnuplot treats "\n\n" as dataset separator (plot ... index n)
    # outputting spec_dry.txt and dry_n to hdf file
    bins = np.empty(self.bins_dry.size - 1)
    for i in range(0, bins.size) :
      prtcls.diag_dry_rng(self.bins_dry[i], self.bins_dry[i+1])
      prtcls.diag_dry_mom(0)
      bins[i] = np.frombuffer(prtcls.outbuf())
    save(self.out_dry, self.bins_dry[0:-1], bins)
    self.dry_n_h5[it_out,:] = bins

    # outputting spec_wet.txt and wet_n to hdf file
    bins = np.empty(self.bins_wet.size - 1)
    for i in range(0, bins.size) :
      prtcls.diag_wet_rng(self.bins_wet[i], self.bins_wet[i+1])
      prtcls.diag_wet_mom(0)
      bins[i] = np.frombuffer(prtcls.outbuf())      
    save(self.out_wet, self.bins_wet[0:-1], bins)
    self.wet_n_h5[it_out,:] = bins
    
    # outputting sounding.txt and sounding.hdf
    self.out_snd.write(u"%g" % (rhod))
    self.out_snd.write(u"\t%g" % (th_d))
    self.out_snd.write(u"\t%g" % (r_v))
    # outputting to hdf 
    self.rho_h5[it_out] = rhod
    self.thd_h5[it_out] = th_d
    self.rv_h5[it_out] =  r_v

    ## cloud water #TODO - hdf
    prtcls.diag_wet_rng(.5e-6, 25e-6)
    for k in range(0,4):
      prtcls.diag_wet_mom(k)
      self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))

    ## chem stuff #TODO import lib
    #prtcls.diag_wet_rng(0,1) # 0 ... 1 m #TODO: consider a select-all option?
    #prtcls.diag_chem(lgrngn.chem_species_t.S_VI)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
    #prtcls.diag_chem(lgrngn.chem_species_t.H)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
    #prtcls.diag_chem(lgrngn.chem_species_t.SO2)
    #self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))

    self.out_snd.write(u"\n")

