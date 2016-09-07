import os
import subprocess as sp
import numpy as np
import kirishima_info as info
import imp
imp.reload(info)



def tephra(results_file=info.results_file,conf_file=info.conf_file,pts_file=info.obs_file,\
		wind_file=info.wind_file,exe_file=info.exe_file):
	args = [exe_file,conf_file,pts_file,wind_file]
	f_out,f_err = open(results_file,'w'),open("stderr.out",'w')
	sp.call(args,stdout=f_out,stderr=f_err)
	f_out.close(),f_err.close()
	return


def pest(case):
    exe = os.path.join("exe","pest")
    args = [exe,case]
    f_out,f_err = open("pest_stdout.out",'w'),open("pest_stderr.out",'w')
    sp.call(args,stdout=f_out,stderr=f_err)
    f_out.close(),f_err.close()
    return    


if __name__ == '__main__':  
    tephra()
    