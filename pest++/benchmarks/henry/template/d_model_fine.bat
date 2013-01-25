@echo off
del ref_coarse\hk.ref
del ref_coarse\vk.ref
del henry_coarse.wel
plproc.exe plproc_coarse.in >nul
par2par.exe par2par_coarse.in >nul
swt_v4x64.exe henry_coarse.nam_swt >nul
mod2obs.exe <mod2obs_head_coarse.in >nul
mod2obs.exe <mod2obs_conc_coarse.in >nul
get_dist_pred.exe <pred_dist_half_henry_coarse.in >nul
get_dist_pred.exe <pred_dist_ten_henry_coarse.in >nul
get_dist_pred.exe <pred_dist_one_henry_coarse.in >nul