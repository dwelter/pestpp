"""10par_xsec tests:
0) "standard user mode" - draw reals from par-bounds prior and obs noise from weights
1) start with existing par csv and obs csv - using empirical parcov and obscov
1a) start with existing par csv and obs csv - using parcov file
2) start with existing par csv and drawing obs en from weights 
3) restart with full simulated obs en
3a) restart with failed runs in simulated obs en
3b) restart with failed runs and bad phi runs in simulated obs en
4) reg_factor = 0.5 test
5) full solution test with standard draw mode
5a) full solution test with empirical parcov

freyberg tests:
0)"standard user mode" - draw reals from par-bounds prior and obs noise from weights
1) draw par en from full parcov supplied in file
2) full solution with empirical parcov - supplied par csv, obs csv and restart csv

synth tests:
0) restart and upgrade 1.1M par problem
"""

