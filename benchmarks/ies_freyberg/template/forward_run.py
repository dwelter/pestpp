import os
import shutil
import pandas as pd
import numpy as np
import pyemu
import flopy
wel_files = ['WEL_0001.dat','WEL_0002.dat']
names=['l','r','c','flux']
wel_fmt = {'l':lambda x: '{0:10.0f}'.format(x)}
wel_fmt['r'] = wel_fmt['l']
wel_fmt['c'] = wel_fmt['l']
wel_fmt['flux'] = lambda x: '{0:15.6E}'.format(x)
for wel_file in wel_files:
    df_w = pd.read_csv(wel_file+'.bak',header=None,names=names,delim_whitespace=True)
    df_t = pd.read_csv(wel_file+'.temp',header=None,names=names,delim_whitespace=True)
    df_t.loc[:,'flux'] *= df_w.flux
    with open(wel_file,'w') as f:
        f.write('        '+df_t.loc[:,names].to_string(index=None,header=None,formatters=wel_fmt)+'\n')
shutil.copy2('WEL_0002.dat','WEL_0003.dat')
pyemu.helpers.run('mfnwt freyberg.nam >_mfnwt.stdout')
pyemu.helpers.run('mp6 freyberg.mpsim >_mp6.stdout')
pyemu.gw_utils.apply_mflist_budget_obs('freyberg.list')
hds = flopy.utils.HeadFile('freyberg.hds')
f = open('freyberg.hds.dat','wb')
for data in hds.get_alldata():
    data = data.flatten()
    np.savetxt(f,data,fmt='%15.6E')
endpoint_file = 'freyberg.mpenpt'
lines = open(endpoint_file, 'r').readlines()
items = lines[-1].strip().split()
travel_time = float(items[4]) - float(items[3])
with open('freyberg.travel', 'w') as ofp:
    ofp.write('travetime {0:15.6e}{1}'.format(travel_time, '\n'))
pyemu.gw_utils.modflow_read_hydmod_file('freyberg.hyd.bin')
