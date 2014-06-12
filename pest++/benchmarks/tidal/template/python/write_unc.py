import os
import numpy as np
import mat_handler as mhand

nrow,ncol = 50,25
delxy = 20.0
row_y = np.arange(0.5*delxy,(nrow+1)*delxy-(0.5*delxy),delxy)
col_x = np.arange(0.5*delxy,(ncol+1)*delxy-(0.5*delxy),delxy)
print len(row_y),len(col_x)
#--write a pplocs file
f = open(os.path.join("misc","pp_locs.dat"),'w')
for i,y in enumerate(row_y):
    for j,x in enumerate(col_x):
        #print i,j
        name = 'p{0:02d}{1:02d}'.format(i+1,j+1)
        f.write("{0:20s} {1:15.6E} {2:15.6E} 1  1.0\n".\
            format(name,x,y))
f.close()

#--write a structure file
a = delxy * 3
aniso = 1.0/100.0
contrib = 1.0
nugget = 0.01
grid_file = os.path.join("misc","xsec_grid.spc")
struct_file = os.path.join("misc","structure.dat")

f = open(struct_file,'w')
f.write("STRUCTURE struct1\n")
f.write(" NUGGET {0:15.6E}\n".format(nugget))
f.write(" TRANSFORM log\n")
f.write(" NUMVARIOGRAM 1\n")
f.write(" VARIOGRAM var1 {0:15.6E}\n".format(contrib))
f.write("END STRUCTURE\n")
f.write("VARIOGRAM var1\n")
f.write(" VARTYPE 2\n")
f.write(" BEARING 0.0\n")
f.write(" A {0:15.6E}\n".format(a))
f.write(" ANISOTROPY {0:15.6E}\n".format(aniso))
f.write("END VARIOGRAM\n")
f.close()
np.savetxt(os.path.join("misc","zone.ref"),np.ones((nrow,ncol)),fmt="%2.0d")
  
args = ["misc\pp_locs.dat","0.0","misc\structure.dat","struct1","unc\cov.dat",'']
f = open(os.path.join("misc","ppcov.in"),'w')
f.write('\n'.join(args)+'\n')
f.close() 

os.system("exe\ppcov.exe <misc\ppcov.in")        
f = open(os.path.join("misc","par.info"),'r')
par_names = {}
for line in f:
    pname = line.strip().split()[0]
    group = pname.split('_')[0]
    if group not in par_names.keys():
        par_names[group] = []
    par_names[group].append(pname)        
cov = mhand.uncert([]) 
var_mult = {"k":1.0,"s":0.1,"sy":0.1}
cov.from_ascii(os.path.join("unc","cov.dat"))
f_unc = open(os.path.join("unc","unc.dat"),'w')
for group,names in par_names.iteritems():
    #print names
    cov.row_names = names
    cov.col_names = names
    cov.to_ascii(os.path.join("unc",group+".dat"))
    f_unc.write("START COVARIANCE_MATRIX\n")
    f_unc.write("FILE "+os.path.join("unc",group+".dat")+'\n')
    f_unc.write("variance_multiplier {0:15.6E}\n".\
        format(var_mult[group[:-2]]))
    f_unc.write("END COVARIANCE_MATRIX\n\n")
f_unc.close()    
                   
        

