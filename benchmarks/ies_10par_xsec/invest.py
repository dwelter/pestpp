import os
import numpy as np

prefixes = ["prior_par_diff", "am_u", "am_s_inv", "obs_diff", "par_diff", "scaled_par_resid", "x4", "x5", "x6", "x7",
            "ivec",".ut","s2","upgrade_1"]

pyemu_dir = "es_pyemu"
pestpp_dir = "master"
pyemu_files = os.listdir(pyemu_dir)
pestpp_files = os.listdir(pestpp_dir)

for prefix in prefixes:
    pyfiles = [f for f in pyemu_files if prefix+'.' in f.lower()]
    #pyarrs = [np.loadtxt(os.path.join(pyemu_dir,f)).flatten() for f in pyfiles]
    for pyf in pyfiles:
        if not pyf in pestpp_files:
            print("pestpp file missing ",pyf)
            continue
        try:
            pyarr = np.loadtxt(os.path.join(pyemu_dir,pyf))
            pparr = np.loadtxt(os.path.join(pestpp_dir,pyf))
        except Exception as e:
            print("error loading arrs",pyf,e)
        if pyarr.shape != pparr.shape:
            print("shape mismatch",pyf,pyarr.shape,pparr.shape)
        try:
            diff = np.abs(pyarr - pparr)
            print(pyf, diff.max(),diff.sum())
        except Exception as e:
            print("error comparing", pyf, e)
