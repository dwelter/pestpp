import os
import numpy as np

nrow,ncol = 80,50
delxy = 10

odtype = np.dtype([('name','a10'),('i','i4'),('j','i4')])
obs = np.genfromtxt('misc\\obs_indices.dat',odtype)  

#--write the coords file
f = open(os.path.join("misc","obs.crd"),'w')

for name,i,j in obs:
    y = (nrow*delxy) - (i * delxy + (0.5 * delxy))
    x = j * delxy + (0.5 * delxy)
    f.write("{0:20s} {1:15.6E} {2:15.6E} 1\n".format(name,x,y))
f.close()    


           


