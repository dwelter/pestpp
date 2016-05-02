import numpy as np

f = open("kirishima_samples.xyz2.utm",'r')
f_ins = open("stage_A.ins",'w')
f_obs = open("pest_obs.dat",'w')
f_ins.write("pif ~\n")
for i,line in enumerate(f):
	val = float(line.strip().split()[-1])
	obs_name = "imas_{0:03d}".format(i+1)
	if i == 0:
		f_ins.write("l2 ["+obs_name+"]67:87\n")
	else:
		f_ins.write("l1 ["+obs_name+"]67:87\n")
	f_obs.write(obs_name+' {0:15.6E}  {1:15.6E} imas\n'.format(val,1.0/val))
f_ins.close()
f_obs.close()	