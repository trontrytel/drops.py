import os
import h5py
import numpy as np
import libcloudphxx as libcl

class output_lgr:

  def __init__(self, outdir, out_time, mom_diag=range(4), chem_sp = ["S_VI", "H", "SO2"],
               cloud_rng = (.5e-6, 25e-6),
               cloud_nbins = 49,
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
    self.out_time = out_time
    self.mom_diag = mom_diag 
    self.chem_sp = chem_sp
    self.cloud_rng = cloud_rng
    self.bins_dry = bins_dry
    self.bins_wet = bins_wet
    self.hdf_spec  = h5py.File(self.outdir_hdf + "spec_drywet.hdf", mode='w')
    self.hdf_sound = h5py.File(self.outdir_hdf + "sounding.hdf", mode='w')
    self.out_snd = open(outdir + "/sounding.txt", mode='w')
    self.out_dry = open(outdir + "/spec_dry.txt", mode='w')
    self.out_wet = open(outdir + "/spec_wet.txt", mode='w')
    self.out_cld = open(outdir + "/spec_cld.txt", mode='w')

    # previous-timestep saturation (RH)
    self.last = {'RH':0, 'cld_mom':{}, 'act_mom':{}}
    self.RH_max = None

    # defining cld bin locations from cloud_rng and cloud_nbins (linear!)
    self.bins_cld = np.linspace(cloud_rng[0], cloud_rng[1], cloud_nbins+1, 
                                endpoint=True)

    # creating data sets that will be used as scales in spec_drywet.hdf               
    self.hdf_spec["time"] = out_time
    self.hdf_spec["bins_dry"] = self.bins_dry[:-1]
    self.hdf_spec["bins_wet"] = self.bins_wet[:-1]
    self.hdf_spec["bins_cld"] = self.bins_cld[:-1]
  
    # creating hdf varibales
    variables_spec = ['dry_number', 'wet_number', "cld_number"]
    # concentration per size of dry particles
    self.dry_number = self.hdf_spec.create_dataset("dry_number",
                           (out_time.size, self.bins_dry.size - 1), dtype='f')
    # concentration per size of wet particles
    self.wet_number = self.hdf_spec.create_dataset("wet_number",
                                 (out_time.size, self.bins_wet.size - 1), dtype='f')
    # concentratio per size of wet particles
    self.cld_number = self.hdf_spec.create_dataset("cld_number",
                                 (out_time.size, self.bins_cld.size - 1), dtype='f')

    # attaching scales, units etc. to the hdf variables
    self.dry_number.dims.create_scale(self.hdf_spec["bins_dry"], "dry_radius")
    self.dry_number.dims[1].attach_scale(self.hdf_spec["bins_dry"])
    self.wet_number.dims.create_scale(self.hdf_spec["bins_wet"], "wet_radius")
    self.wet_number.dims[1].attach_scale(self.hdf_spec["bins_wet"])
    self.cld_number.dims.create_scale(self.hdf_spec["bins_cld"], "cld_radius")
    self.cld_number.dims[1].attach_scale(self.hdf_spec["bins_cld"])

    for var in variables_spec:
      getattr(self,var).dims.create_scale(self.hdf_spec["time"], "time")
      getattr(self,var).dims[0].attach_scale(self.hdf_spec["time"])
      getattr(self,var).dims[0].label = 's'
      getattr(self,var).dims[1].label = 'm'
      getattr(self,var).attrs["Units"] = "1/kg"


    # creating data sets that will be used as scales in sounding.hdf
    self.hdf_sound["time"] = out_time


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
      #print "var sound", var
      #pdb.set_trace()
      setattr(self, var, self.hdf_sound.create_dataset(var, (out_time.size,), dtype='f')) 
      getattr(self, var).dims.create_scale(self.hdf_sound["time"], "time")
      getattr(self, var).dims[0].attach_scale(self.hdf_sound["time"])
      getattr(self, var).dims[0].label = 's'
      getattr(self, var).attrs["Units"] = units_sound[var]


    # description of the txt files
    self.out_snd.write(u"#rhod [kg/m3]\tth_d [K] (theta dry!)\tr_v [kg/kg] (mixing ratio)\tM0 [TODO]\tM1 [TODO]\tM2 [TODO]\tM3 [TODO]\n")
    self.out_dry.write(u"#r_d [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
    self.out_wet.write(u"#r_w [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
    self.out_cld.write(u"#r_w [m] (left bin edge)\tn [kg-1] (per mass of dry air)\n")
  
  
    # placing a quick-look gnuplot file in the output directory
    import shutil
    shutil.copyfile(os.path.dirname(__file__) + '/quicklook.gpi', outdir + '/quicklook.gpi')

  def diag(self, prtcls, 
    # values before adjustment by microphysics
    rhod_in, th_d_in, r_v_in, 
    # values after adjustment by microphysics
    rhod, th_d, r_v, 
    itime, 
    save = True, stats = None
  ):
    def S(rhod, th_d, r_v):
      T  = libcl.common.T(th_d[0], rhod[0])
      return (rhod * r_v) * libcl.common.R_v * T / libcl.common.p_vs(T)

    RH_in = S(rhod_in, th_d_in, r_v_in)
    RH = S(rhod, th_d, r_v)

    if self.RH_max is None and RH < self.last['RH']:
      #print self.last['RH_in'], self.last['RH'], RH_in, RH
      stats['S_max_RH'] = max(RH_in, self.last['RH_in'])
      self.RH_max = stats['S_max_RH']
      stats['S_max_A0'] = self.last['act_mom'][0]
      stats['S_max_M0'] = self.last['cld_mom'][0]
      stats['S_max_M1'] = self.last['cld_mom'][1]
      stats['S_max_M2'] = self.last['cld_mom'][2]
      stats['S_max_M3'] = self.last['cld_mom'][3]
      stats['S_max_rhod'] = self.last['rhod']

    self.last['RH'] = RH
    self.last['RH_in'] = RH_in
    self.last['rhod'] = rhod
      
    ## cloud water 
    prtcls.diag_wet_rng(self.cloud_rng[0], self.cloud_rng[1])
    for k in self.mom_diag: 
      prtcls.diag_wet_mom(k)
      self.last['cld_mom'][k] = np.frombuffer(prtcls.outbuf())[0]

    ## activated droplets
    prtcls.diag_rw_ge_rc()
    prtcls.diag_wet_mom(0)
    self.last['act_mom'][0] = np.frombuffer(prtcls.outbuf())[0]

    ## all what's below happens only every outfreq timesteps
    if not save:
      return

    # index of itime in the array self.time
    it_out = self.out_time.searchsorted(itime)

    def bins_save(type):
      def save(out, xx, yy):
	for x, y in zip(xx, yy):
	  out.write(u"%g\t%g\n" % (x, y))
	out.write(u"\n\n") # gnuplot treats "\n\n" as dataset separator (plot ... index n)

      bins = np.empty(getattr(self, 'bins_' + type).size - 1)
      for i in range(bins.size) :
        if (type in ['cld', 'wet']):
          prtcls.diag_wet_rng(getattr(self, 'bins_' + type)[i], getattr(self, 'bins_' + type)[i+1])
          prtcls.diag_wet_mom(0)
        elif (type in ['dry']):
          prtcls.diag_dry_rng(getattr(self, 'bins_' + type)[i], getattr(self, 'bins_' + type)[i+1])
          prtcls.diag_dry_mom(0)
        bins[i] = np.frombuffer(prtcls.outbuf())
      save(getattr(self, 'out_' + type), getattr(self, 'bins_' + type)[0:-1], bins)
      getattr(self, type + '_number')[it_out,:] = bins

    for type in ["cld", "wet", "dry"]:
      bins_save(type)
    
    # outputting sounding.txt and sounding.hdf
    self.out_snd.write(u"%.9g" % (rhod))
    self.out_snd.write(u"\t%.9g" % (th_d))
    self.out_snd.write(u"\t%.9g" % (r_v))
    # outputting to hdf 
    self.rhod[it_out] = rhod
    self.thd[it_out] = th_d
    self.rv[it_out] =  r_v

    ## cloud water 
    for k in self.mom_diag: 
      self.out_snd.write(u"\t%.9g" % (self.last['cld_mom'][k]))
      getattr(self,'mom_' + str(k))[it_out] = self.last['cld_mom'][k]

    ## chem stuff 
    prtcls.diag_wet_rng(0,1) # 0 ... 1 m #TODO: consider a select-all option?
    for sp in self.chem_sp:
      prtcls.diag_chem(getattr(libcl.lgrngn.chem_species_t, sp)) 
      getattr(self, "conc_" + sp)[it_out] =  np.frombuffer(prtcls.outbuf())
 
    self.out_snd.write(u"\n")
