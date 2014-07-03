import os
import shutil
import socket
import subprocess as sp

BENCHMARK_DIRS = ["3pg","10par_xsec","box","henry","kirishima","stor"]#,"tidal"]
#BENCHMARK_DIRS = ["henry"]
NSLAVES = 4
BEOPEST_EXE = os.path.join("exe","beopest.exe")
PSTNAME = "pest.pst"
HOSTNAME = socket.gethostname()
PESTPP_EXE = os.path.join("..","x64","release","pest++.exe")
def beopest():
	for bdir in BENCHMARK_DIRS:
		template_dir = os.path.join(bdir,"template")
		beopest_dir = os.path.join(bdir,"beopest")
		if os.path.exists(beopest_dir):
			try:
				shutil.rmtree(beopest_dir)
			except Exception as e:
				print "unable to remove existing beopest dir: ",beopest_dir
				raise Exception(e)
		os.mkdir(beopest_dir)
		
		slave_dir = os.path.join(beopest_dir,"slaves")
		if os.path.exists(slave_dir):
			try:
				shutil.rmtree(master_dir)
			except Exception as e:
				print "unable to remove existing beopest slave dir: ",slave_dir
				raise Exception(e)
		os.mkdir(slave_dir)	

		master_dir = os.path.join(beopest_dir,"master")		
		if os.path.exists(master_dir):
			try:
				shutil.rmtree(master_dir)
			except Exception as e:
				print "unable to remove existing beopest master dir: ",master_dir
				raise Exception(e)
		shutil.copytree(template_dir,master_dir)

		#--start the master
		os.chdir(master_dir)
		sp.Popen([BEOPEST_EXE,PSTNAME,"/h",":4004"],shell=False)
		os.chdir(os.path.join("..","..",".."))
		

		#--start the slaves
		os.chdir(slave_dir)

		args = [os.path.join("..","..","..","hpc_client_util.exe"),"-cmdExec:"+BEOPEST_EXE,\
		    "-cmdArgs:\" "+PSTNAME+" /H "+HOSTNAME+":4004 \"",\
			"-src:"+os.path.join("..","..","template"),"-n:"+str(NSLAVES)]					
		print args
		sp.call(args,shell=False)
		os.chdir(os.path.join("..","..",".."))


def pestpp():
	for bdir in BENCHMARK_DIRS:
		template_dir = os.path.join(bdir,"template")
		pestpp_dir = os.path.join(bdir,"pestpp")
		if os.path.exists(pestpp_dir):
			try:
				shutil.rmtree(pestpp_dir)
			except Exception as e:
				print "unable to remove existing pestpp dir: ",pestpp_dir
				raise Exception(e)
		os.mkdir(pestpp_dir)
		
		slave_dir = os.path.join(pestpp_dir,"slaves")
		if os.path.exists(slave_dir):
			try:
				shutil.rmtree(master_dir)
			except Exception as e:
				print "unable to remove existing pestpp slave dir: ",slave_dir
				raise Exception(e)
		os.mkdir(slave_dir)	

		master_dir = os.path.join(pestpp_dir,"master")		
		if os.path.exists(master_dir):
			try:
				shutil.rmtree(master_dir)
			except Exception as e:
				print "unable to remove existing pestpp master dir: ",master_dir
				raise Exception(e)
		
		#shutil.copy2(PESTPP_EXE,os.path.join(template_dir,"pest++.exe"))
		shutil.copytree(template_dir,master_dir)

		#--start the master
		os.chdir(master_dir)
		pestpp_exe = os.path.join("..","..","..",PESTPP_EXE)
		print pestpp_exe
		sp.Popen([pestpp_exe,PSTNAME,"/h",":4004"],shell=False)
		os.chdir(os.path.join("..","..",".."))
	
		#--start the slaves
		os.chdir(slave_dir)
		pestpp_exe = os.path.join("..","..","..","..",PESTPP_EXE)
		args = [os.path.join("..","..","..","hpc_client_util.exe"),"-cmdExec:"+pestpp_exe,\
		    "-cmdArgs: \' "+PSTNAME+" /H "+HOSTNAME+":4004 \'",\
			"-src:"+os.path.join("..","..","template"),"-n:"+str(NSLAVES)]					
		print args
		sp.call(args,shell=False)
		os.chdir(os.path.join("..","..",".."))
	

def plot_phi_vs_iter():
	plt_dir = os.path.join("plots","phi")
	if os.path.exists(plt_dir):
		shutil.rmtree(plt_dir)
	os.makedirs(plt_dir)
	for bdir in BENCHMARK_DIRS:
		print bdir

		beopest_dir = os.path.join(bdir,"beopest","master")
		beopest_rec = os.path.join(beopest_dir,"pest.rec")
		assert os.path.exists(beopest_rec), "missing beopest rec file: "+str(beopest_rec)

		pestpp_dir = os.path.join(bdir,"pestpp","master")
		pestpp_rec = os.path.join(pestpp_dir,"pest.rec")
		assert os.path.exists(pestpp_rec), "missing pestpp rec file: "+str(pestpp_rec)

		beopest_phi = get_phi_beopest(beopest_rec)
		pestpp_phi = get_phi_pestpp(pestpp_rec)
		print beopest_phi[:,2]
		print pestpp_phi[:,2]

		fig = pylab.figure()
		ax = pylab.subplot(111)
		axt = pylab.twinx()
		l1 = ax.plot(beopest_phi[:,0],beopest_phi[:,2],'b--',marker='.',label="phi - beopest")
		l2 = axt.plot(beopest_phi[:,0],beopest_phi[:,1],'g--',label="model calls - beopest")
		l3 = ax.plot(pestpp_phi[:,0],pestpp_phi[:,2],'b-',marker='.',label="phi - pestpp")
		l4 = axt.plot(pestpp_phi[:,0],pestpp_phi[:,1],'g-',label="model calls - pestpp")
		ax.grid()
		axt.set_ylabel("model calls")
		ax.set_ylabel("phi")
		ax.set_xlabel("iteration number")
		lines = l1+l2+l3+l4
		labs = [l.get_label() for l in lines]
		ax.legend(lines,labs,loc=1,ncol=2)
		pylab.savefig(os.path.join(plt_dir,bdir+".png"))

def get_phi_beopest(rec_file):
	beo_key = "OPTIMISATION ITERATION NO."
	phi_rec = []
	f = open(rec_file,'r')
	while True:
		line = f.readline()
		if line == '':
			break
		if beo_key in line:
			raw = line.strip().split()
			inum = int(raw[-1])
			line = f.readline()
			raw = line.strip().split()
			model_calls = int(raw[-1])
			line = f.readline()
			if "regul" in line:
				line = f.readline()
			raw = line.strip().split()
			phi = float(raw[-1])
			phi_rec.append([inum,model_calls,phi])
	f.close()
	return np.array(phi_rec)


def get_phi_pestpp(rec_file):
	beo_key = "OPTIMISATION ITERATION NUMBER"
	phi_rec = []
	f = open(rec_file,'r')
	while True:
		line = f.readline()
		if line == '':
			break
		if beo_key in line:
			raw = line.strip().split()
			inum = int(raw[-1])
			for _ in range(5):
				line = f.readline()
			raw = line.strip().split()
			model_calls = int(raw[-1])
			line = f.readline()
			line = f.readline()
			raw = line.strip().split()
			phi = float(raw[-1])
			phi_rec.append([inum,model_calls,phi])
	f.close()
	return np.array(phi_rec)





if __name__ == "__main__":
	beopest()
	#pestpp()


