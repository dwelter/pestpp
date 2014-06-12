import os
import numpy as np

plot_dir = "plots"
if not os.path.exists(plot_dir):
	os.mkdir(plot_dir)

vent_x,vent_y = 110616.0,3538377.0
map_xlim = [vent_x-10000, vent_x+50000]
map_ylim = [vent_y-50000,vent_y+10000]
grid_step = 250

obs_file = os.path.join('background','kirishima_samples.xyz2.utm')
exe_file = os.path.join("exe","tephra2")
conf_file = os.path.join("conf","stage_A.conf")
wind_file = os.path.join("wind","stage_A.wind")
results_file = "stage_A.results"
grid_file=os.path.join("background","kirishima_grid.xy")
grid_results_file=os.path.join("background","kirishima_grid.xy.results")


obs_xys = np.loadtxt(obs_file, usecols=[0,1,3])
min_mass = 1.0


def write_grid():
    #--write the grid xy file
    xs = np.arange(map_xlim[0],map_xlim[1],grid_step)
    ys = np.arange(map_ylim[0],map_ylim[1],grid_step)
    f = open(grid_file,'w')
    for x in xs:
        for y in ys:
            f.write("{0:15.6E}  {1:15.6E}  0.0  \n".format(x,y))
    f.close()


def  resize(x,y,data):
    unique_x = np.sort(np.unique(x))
    unique_y = np.sort(np.unique(y))
    nrow,ncol = unique_y.shape[0],unique_x.shape[0]
    i_idx,j_idx = {},{}
    for j,ux in enumerate(unique_x):
        j_idx[ux] = j
    for i,uy in enumerate(unique_y):
        i_idx[uy] = i
    arr = np.zeros((nrow,ncol))-999
    for x,y,val in zip(x,y,data):        
        arr[i_idx[y],j_idx[x]] = val
    return arr


def load_rei(case):
    #--get the number of observation from the pst - use to cut out regul obs
    f = open(case+".pst")
    for iline in range(4):
        line = f.readline()
    f.close()    
    nobs = int(line.strip().split()[1])
    #--load the residuals from each iteration
    if 'reg' in case:
        skip = 7
    else:
        skip = 5
    rei_dtype = np.dtype([('name', 'a12'), ('measured', np.float64),\
                     ('modeled', np.float64),('weight',np.float64)])
    files = os.listdir(os.path.join('.'))
    rei_arrays = {}
    for f in files:
        raw = f.split('.')
        if len(raw) == 3 and raw[0] == case and raw[1] == 'rei':
            iter = int(raw[-1])
            rei = np.genfromtxt(f, skiprows=skip,\
             usecols=[0, 2, 3, 5], dtype=rei_dtype)
            rei_arrays[iter] = rei[:nobs]
    if len(rei_arrays.keys()) == 0:
        if case+".rei" in files:
            rei = np.genfromtxt(case+".rei", skiprows=skip,\
            usecols=[0, 2, 3, 5], dtype=rei_dtype)
            rei_arrays[0] = rei[:nobs]            
    return rei_arrays    