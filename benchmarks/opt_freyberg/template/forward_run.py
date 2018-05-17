import os
import numpy as np
import pandas as pd
import flopy
import pyemu
try:
   os.remove('freyberg.list')
except Exception as e:
   print('error removing tmp file:freyberg.list')
try:
   os.remove('freyberg.hds')
except Exception as e:
   print('error removing tmp file:freyberg.hds')
try:
   os.remove('mt3d_link.ftl')
except Exception as e:
   print('error removing tmp file:mt3d_link.ftl')
try:
   os.remove('freyberg.mt3d')
except Exception as e:
   print('error removing tmp file:freyberg.mt3d')
try:
   os.remove('MT3D001.ucn')
except Exception as e:
   print('error removing tmp file:MT3D001.ucn')
pyemu.helpers.apply_array_pars()

pyemu.gw_utils.apply_sfr_seg_parameters()
pyemu.helpers.run('mfnwt freyberg.truth.nam 1>freyberg.truth.nam.stdout 2>freyberg.truth.nam.stderr')
pyemu.helpers.run('mt3dusgs freyberg.mt3d.nam')

pyemu.gw_utils.apply_mflist_budget_obs('freyberg.list',flx_filename='flux.dat',vol_filename='vol.dat',start_datetime='01-01-1970')
pyemu.gw_utils.apply_hds_obs('freyberg.hds')
pyemu.gw_utils.apply_mtlist_budget_obs('freyberg.mt3d.list')
pyemu.gw_utils.apply_sft_obs()
pyemu.gw_utils.apply_hds_obs('MT3D001.ucn')

