import os
import h5py
import numpy as np
import libcloudphxx as libcl

class output_lgr:

  def __init__(self, outdir, time_out, mom_diag=range(4), chem_sp = ["S_VI", "H", "SO2"],
               cloud_rng = (.5e-6, 25e-6),
               bins_dry = 1e-6 * pow(10, -3 + np.arange(40) * .1),
               bins_wet = 1e-6 * pow(10, -3 + np.arange(25) * .2)):
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
    self.mom_diag = mom_diag 
    self.chem_sp = chem_sp
    self.cloud_rng = cloud_rng
    self.bins_dry = bins_dry
    self.bins_wet = bins_wet
    self.hdf_spec = h5py.File(self.outdir_hdf + "spec_drywet.hdf", mode='w')
    self.hdf_sound = h5py.File(self.outdir_hdf + "sounding.hdf", mode='w')
    self.out_snd = open(outdir + "/sounding.txt", mode='w')
    self.out_dry = open(outdir + "/spec_dry.txt", mode='w')
    self.out_wet = open(outdir + "/spec_wet.txt", mode='w')


    # creating data sets that will be used as scales in spec_drywet.hdf               
    self.hdf_spec["time"] = self.time
    self.hdf_spec["bins_dry"] = self.bins_dry[:-1]
    self.hdf_spec["bins_wet"] = self.bins_wet[:-1]
  
  # creating hdf varibales
    variables_spec = ['dry_number', 'wet_number']
    # concentratio per size of dry particles
    self.dry_number = self.hdf_spec.create_dataset("dry_number",
                                 (self.time.size, self.bins_dry.size - 1), dtype='f')
    # concentratio per size of wet particles
    self.wet_number = self.hdf_spec.create_dataset("wet_number",
                                 (self.time.size, self.bins_wet.size - 1), dtype='f')

    # attaching scales, units etc. to the hdf variables
    self.dry_number.dims.create_scale(self.hdf_spec["bins_dry"], "dry_radius")
    self.dry_number.dims[1].attach_scale(self.hdf_spec["bins_dry"])
    self.wet_number.dims.create_scale(self.hdf_spec["bins_wet"], "wet_radius")
    self.wet_number.dims[1].attach_scale(self.hdf_spec["bins_wet"])
    for var in variables_spec:
      getattr(self,var).dims.create_scale(self.hdf_spec["time"], "time")
      getattr(self,var).dims[0].attach_scale(self.hdf_spec["time"])
      getattr(self,var).dims[0].label = 's'
      getattr(self,var).dims[1].label = 'm'
      getattr(self,var).attrs["Units"] = "1/kg"


    # creating data sets that will be used as scales in sounding.hdf
    self.hdf_sound["time"] = self.time

    # creating hdf varibales and attaching scales, units etc.  
    # sounding: dry density, dry pot. temp., water vapour mix. rat.
    variables_sound = ["rhod", "thd", "rv"]
    units_sound = {"rhod":"kg/m^3", "thd":"K", "rv":"kg/kg"} 
    for i in self.mom_diag:
      variables_sound.append("mom_" + str(i))
      units_sound["mom_" + str(i)] = "m^" + str(i) #TODO is it ok??
    for sp in self.chem_sp:
      variables_sound.append("conc_" + str(sp)) #TODO is it indeed concentration?
      units_sound["conc_" + str(sp)] = "TODO" 

    for var in variables_sound:
      print "var sound", var
      #pdb.set_trace()
      setattr(self, var, self.hdf_sound.create_dataset(var, (self.time.size,), dtype='f')) 
      getattr(self, var).dims.create_scale(self.hdf_sound["time"], "time")
      getattr(self, var).dims[0].attach_scale(self.hdf_sound["time"])
      getattr(self, var).dims[0].label = 's'
      getattr(self, var).attrs["Units"] = units_sound[var]


    # description of the txt files #TODO removing chem?
    self.out_snd.write(u"#rhod [kg/m3]\tth_d [K] (theta dry!)\tr_v [kg/kg] (mixing ratio)\tM0 [TODO]\tM1 [TODO]\tM2 [TODO]\tM3 [TODO]\tS_VI [kg/kg]\tH [kg/kg]\tSO2 [kg/kg]\n")
    self.out_dry.write(u"#r_d [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
    self.out_wet.write(u"#r_w [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")


  def diag(self, prtcls, rhod, th_d, r_v, itime):
    # index of itime in the array self.time
    it_out = self.time.searchsorted(itime)

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
    self.dry_number[it_out,:] = bins

    # outputting spec_wet.txt and wet_n to hdf file
    bins = np.empty(self.bins_wet.size - 1)
    for i in range(0, bins.size) :
      prtcls.diag_wet_rng(self.bins_wet[i], self.bins_wet[i+1])
      prtcls.diag_wet_mom(0)
      bins[i] = np.frombuffer(prtcls.outbuf())      
    save(self.out_wet, self.bins_wet[0:-1], bins)
    self.wet_number[it_out,:] = bins
    
    # outputting sounding.txt and sounding.hdf
    self.out_snd.write(u"%g" % (rhod))
    self.out_snd.write(u"\t%g" % (th_d))
    self.out_snd.write(u"\t%g" % (r_v))
    # outputting to hdf 
    self.rhod[it_out] = rhod
    self.thd[it_out] = th_d
    self.rv[it_out] =  r_v

    ## cloud water 
    prtcls.diag_wet_rng(self.cloud_rng[0], self.cloud_rng[1])
    for k in self.mom_diag: 
      prtcls.diag_wet_mom(k)
      self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
      getattr(self,'mom_' + str(k))[it_out] = np.frombuffer(prtcls.outbuf())

    ## chem stuff 
    prtcls.diag_wet_rng(0,1) # 0 ... 1 m #TODO: consider a select-all option?
    for sp in self.chem_sp:
      prtcls.diag_chem(getattr(libcl.lgrngn.chem_species_t, sp)) 
      self.out_snd.write(u"\t%g" % (np.frombuffer(prtcls.outbuf())))
      getattr(self, "conc_" + sp)[it_out] =  np.frombuffer(prtcls.outbuf())
 
    self.out_snd.write(u"\n")

