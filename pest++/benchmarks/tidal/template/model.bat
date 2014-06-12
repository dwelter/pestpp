
cd model
MODFLOW-NWT_64.exe Tidal_MultiLayer.nam
cd ..

exe\mod2smp.exe <misc\mod2smp_hds.in
exe\tsproc.exe <misc\tsproc_model_run.in

