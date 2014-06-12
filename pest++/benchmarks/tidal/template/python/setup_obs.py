import os
import numpy as np

nrow,ncol = 50,25
obs_col = 12
obs_lays = [1,3]
obs_start_row = 48
nobs = 20
nobs_row = 1
delxy = 20

obs_x = obs_col * delxy - (0.5*delxy)
obs_dely = nobs_row * delxy
obs_start_y = (nrow - obs_start_row) * delxy - (0.5*delxy)
obs_end_y = obs_start_y + (obs_dely * nobs)
obs_y = np.arange(obs_start_y,obs_end_y,obs_dely)


#--write the grid file
f = open(os.path.join("misc","grid.spc"),'w')
f.write("{0:d} {1:d}\n".format(nrow,ncol))
f.write("{0:10.3f} {1:10.3f} 0.0\n".format(0.0,nrow*delxy))
f.write("{0:d}*{1:d}\n".format(ncol,delxy))
f.write("{0:d}*{1:d}\n".format(nrow,delxy))
f.close()

#--write the coords file
f = open(os.path.join("misc","obs.crd"),'w')

for lay in obs_lays:
    for iy,y in enumerate(obs_y):
        name = "o{0:1d}{1:1d}".format(lay,iy)
        f.write("{0:20s} {1:15.6e} {2:15.6e} {3:d}\n"\
            .format(name,obs_x,y,lay))
f.close()


           


