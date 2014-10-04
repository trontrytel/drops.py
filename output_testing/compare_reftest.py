import numpy as np
import h5py
import os, sys, getopt
import pdb

def comparing_hdf(dir_new, dir_ref, filetype='hdf'):
    for filename in os.listdir(dir_ref+"/hdf_output/"):
        if filename.endswith(filetype):
            print "\n", filename
            f_ref = h5py.File(dir_ref+"/hdf_output/" + filename, "r")
            f_new = h5py.File(dir_new+"/hdf_output/" + filename, "r")
            for key in f_ref.keys():
                if (np.array(f_ref[key]) == np.array(f_new[key])).all():
                    print "\n", "values of " +  key + " OK "
                    if key not in ["time", "bins_dry", "bins_wet", "mom_ord", "chem_sp"]:
                        print "new files: unit of ", key, "-", f_new[key].attrs["Units"]
                        print "new files: dimensions of ", key, ":"
                        for dim in f_new[key].dims:
                            print "name, unit, values", dim[0].name, dim.label, dim[0][:]
                    
                else:
                     print "Problems!! values of " + key + " differ"
            
                


def main(dir_new, dir_ref):
    comparing_hdf(dir_new, dir_ref)


if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "n:r:",
                               ["dir_new=", "dir_ref="])
    arg_dic = {}
    for opt, arg in opts:
        if opt in ("-n", "--dir_new"):
            arg_dic["dir_new"] = arg
        elif opt in ("-r", "--dir_ref"):
            arg_dic["dir_ref"] = arg
       
    if "dir_ref" not in arg_dic:
        print "please specify the reference directory, e.g. --dir_ref test_ref"
    elif "dir_new" not in arg_dic:
        print "please specify the new directory you want to compare, e.g. --dir_new ../test"
    else:
        main(**arg_dic)
