import os
import shutil
import platform
import numpy as np
import pyemu


if "windows" in platform.platform().lower():
    exe_path = os.path.join("..", "..", "..", "exe", "windows", "x64", "Release", "pestpp-opt.exe")
elif "darwin" in platform.platform().lower():
    exe_path = os.path.join("..", "..", "..", "exe", "mac", "pestpp-opt")
else:
    exe_path = os.path.join("..", "..", "..", "exe", "linux", "pestpp-opt")

def std_weights_test():
	d = os.path.join("opt_dewater_chance","test_std_weights2")
	if os.path.exists(d):
		shutil.rmtree(d)
	shutil.copytree(os.path.join("opt_dewater_chance","template"),d)
	pst_file = os.path.join(d,"dewater_pest.base.pst")
	jco_file = os.path.join(d,"dewater_pest.full.jcb")
	pst = pyemu.Pst(pst_file)
	par = pst.parameter_data
	par.loc[par.partrans=="fixed","partrans"] = "log"
	jco = pyemu.Jco.from_binary(jco_file)
	par.loc[par.pargp=="q","partrans"] = "fixed"
	obs = pst.observation_data.loc[pst.nnz_obs_names,:]

	forecast_names = list(obs.loc[obs.obgnme.apply(lambda x: x.startswith("l_") or\
							 x.startswith("less_")),"obsnme"])
	#print(forecast_names)
	pst.observation_data.loc[:,"weight"] = 0.0
	sc = pyemu.Schur(jco=jco,pst=pst,forecasts=forecast_names)
	#print(sc.get_forecast_summary())
	
	fstd = sc.get_forecast_summary().loc[:,"post_var"].apply(np.sqrt)	
	pr_unc_py = fstd.to_dict()
	pst.observation_data.loc[fstd.index,"weight"] = fstd.values

	pst.pestpp_options["opt_risk"] = 0.1
	pst.pestpp_options["opt_std_weights"] = True
	pst.pestpp_options["base_jacobian"] = os.path.split(jco_file)[-1]
	par.loc[par.pargp=="q","partrans"] = "none"
	
	new_pst_file = os.path.join(d,"test.pst")

	pst.write(new_pst_file)
	pyemu.os_utils.run("{0} {1}".format(exe_path,os.path.split(new_pst_file)[-1]),cwd=d)
	pr_unc1 = scrap_rec(new_pst_file.replace(".pst",".rec"))


	print(pr_unc1)
	for fore in forecast_names:
		dif = np.abs(pr_unc_py[fore]-pr_unc1[fore])
		print(fore,pr_unc_py[fore],pr_unc1[fore],dif)
		assert dif < 1.0e-4

	pst.pestpp_options["opt_std_weights"] = False
	pst.write(new_pst_file)
	pyemu.os_utils.run("{0} {1}".format(exe_path,os.path.split(new_pst_file)[-1]),cwd=d)
	pr_unc2 = scrap_rec(new_pst_file.replace(".pst",".rec"))
	print(pr_unc2)
	for fore in forecast_names:
		dif = np.abs(pr_unc_py[fore]-pr_unc2[fore])
		print(fore,pr_unc_py[fore],pr_unc2[fore],dif)
		assert dif < 1.0e-4


def scrap_rec(rec_file):
	unc = {}
	tag = "FOSM-based chance constraint information at start of iteration 1"
	with open(rec_file,'r') as f:
		while True:
			line = f.readline()
			if line == "":
				break
			
			if tag in line:
				f.readline()
				while True:
					line = f.readline()
					if line == "":
						break
					raw = line.strip().split()
					#print(raw)
					try:
						name = raw[0].lower()
						val = float(raw[4])
						unc[name] = val
					except:
						break
				break
	return unc

def basic_run_test():
	pass

if __name__ == "__main__":
	std_weights_test()