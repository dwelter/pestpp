import os
import shutil
import socket
import subprocess as sp

BENCHMARK_DIRS = ["3pg","10par_xsec","box","henry","kirishima","stor","tidal"]
NSLAVES = 6
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
		    "-cmdArgs: \" "+PSTNAME+" /H "+HOSTNAME+":4004 \"",\
			"-src:"+os.path.join("..","..","template"),"-n:"+str(NSLAVES)]					
		print args
		sp.call(args,shell=False)
		os.chdir(os.path.join("..","..",".."))
		break


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
		    "-cmdArgs: \" "+PSTNAME+" /H "+HOSTNAME+":4004 \"",\
			"-src:"+os.path.join("..","..","template"),"-n:"+str(NSLAVES)]					
		print args
		sp.call(args,shell=False)
		os.chdir(os.path.join("..","..",".."))
		break


#def pestpp():
#	if not os.path.exists("slaves"):
#		os.mkdir("slaves")	
#
#	for bdir in BENCHMARK_DIRS:
#		template_dir = os.path.join(bdir,"template")
#		master_dir = os.path.join(bdir,"pestpp_master")
#		if os.path.exists(master_dir):
#			try:
#				shutil.rmtree(master_dir)
#			except Exception as e:
#				print "unable to remove existing pestpp master dir: ",master_dir
#				raise Exception(e)
#		shutil.copytree(template_dir,master_dir)
#		
#		#--start the master
#		os.chdir(master_dir)
#		sp.Popen(["pest++.exe",PSTNAME,"/h",":4004"],shell=False)
#		os.chdir(os.path.join("..",".."))
#
#		#--start the slaves
#		args = [os.path.join("..","hpc_client_util.exe"),"-cmdExec:pest++.exe",\
#		    "-cmdArgs: \" "+PSTNAME+" /H "+HOSTNAME+":4004 \"",\
#			"-src:"+os.path.join("..",template_dir),"-n:"+str(NSLAVES)]
#				
#		os.chdir("slaves")
#		sp.call(args,shell=False)
#		os.chdir(os.path.join(".."))
#		#break





if __name__ == "__main__":
	beopest()
	pestpp()


