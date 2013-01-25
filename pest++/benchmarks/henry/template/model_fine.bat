@echo off
del ref_fine\hk.ref >nul
del ref_fine\vk.ref >nul
del henry_fine.wel >nul
plproc.exe plproc_fine.in >nul
par2par.exe par2par_fine.in >nul
swt_v4x64.exe henry_fine.nam_swt >nul
mod2obs.exe <mod2obs_head_fine.in >nul
mod2obs.exe <mod2obs_conc_fine.in >nul
get_dist_pred.exe <pred_dist_half_henry_fine.in >nul
get_dist_pred.exe <pred_dist_ten_henry_fine.in  >nul
get_dist_pred.exe <pred_dist_one_henry_fine.in  >nul