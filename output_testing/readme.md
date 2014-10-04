Simple python code to compare output to the reference one. 

Results from 2 simulations are included: short one (1 time step) - ``test_nt1_ref``, and longer one - ``test_ref``.
Options used in simulations are in readme files in both directories.

Examples of usage:

    $ python compare_reftest.py --dir_ref test_nt1_ref --dir_new ../test_nt1
  
    $ python compare_reftest.py --dir_ref test_ref --dir_new ../test