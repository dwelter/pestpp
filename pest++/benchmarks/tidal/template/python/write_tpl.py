import os
import numpy as np

nrow,ncol = 50,25
par_dict = {'k':[1,2,3],'s':[2,3],'sy':[1]}
par_values = {'k':1000,'s':1.0e-5,'sy':0.15}    
ref_dir = os.path.join("model","ref")
tpl_dir = "tpl"
if not os.path.exists(tpl_dir):
    os.mkdir(tpl_dir)
base_arr = np.zeros((nrow,ncol))
par_names = {}
for par,layers in par_dict.iteritems():
    par_names[par] = []
    for lay in layers:
        name = par+"_{0:d}.ref".format(lay)
        np.savetxt(os.path.join(ref_dir,name),base_arr+par_values[par],fmt="%15.6E")        
        f_tpl = open(os.path.join(tpl_dir,name+".tpl"),'w')
        f_tpl.write("ptf ~\n")
        for i in range(nrow):
            for k in range(ncol):
                pname = par+"{0:02d}_{1:02d}_{2:02d}".format(lay,i+1,k+1)
                f_tpl.write("~ {0:18s} ~ ".format(pname))
                par_names[par].append(pname)
            f_tpl.write('\n')
        f_tpl.close()        
                        
f_par = open(os.path.join("misc","par.info"),'w')
f_grp = open(os.path.join("misc","grp.info"),'w')
for par,pnames in par_names.iteritems():
    f_grp.write(par+" relative     1.0000E-02   0.000      switch      2.000      parabolic\n")
    for pname in pnames:
        #print par,pname,par_values[par]
        f_par.write("{0:20s}  log factor  {1:15.6E}  1.0e-10 1.0e+10 {2:10s} 1 0 1\n".\
            format(pname,par_values[par],par))  

f_grp.close()
f_par.close()            



