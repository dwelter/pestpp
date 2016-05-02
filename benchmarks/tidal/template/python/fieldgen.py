import os
import shutil
import numpy as np
import pylab

nrow,ncol = 50,25

np.savetxt(os.path.join("misc","zone.ref"),np.ones((nrow,ncol)),fmt="%3d")
real_dir = "reals"
if os.path.exists(real_dir):
	shutil.rmtree(real_dir)
os.mkdir(real_dir)
args = [os.path.join("misc","grid.spc"),'',\
		os.path.join("misc","zone.ref"),'f',\
		os.path.join("misc","structure.dat"),\
		"struct1",'o',"10","10",os.path.join(real_dir,"real"),'f',
		"1000",'']
f = open(os.path.join("misc","fieldgen.in"),'w')
f.write('\n'.join(args)+'\n')
f.close()
os.system("exe\\fieldgen.exe  <misc\\fieldgen.in")

real_files = os.listdir(real_dir)
for real_file in real_files:
	arr = np.log10(np.loadtxt(os.path.join(real_dir,real_file)))
	pylab.imshow(arr,interpolation="none")
	print real_file	
	pylab.show()
