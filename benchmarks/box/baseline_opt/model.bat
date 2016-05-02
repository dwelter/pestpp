@echo off
del model\ref\hk_layer_1.ref
exe\fac2real.exe <misc\fac2real.in >nul

cd model
mf2005.exe syn.nam >nul
cd ..

exe\mod2smp.exe <misc\mod2smp_hds.in >nul
