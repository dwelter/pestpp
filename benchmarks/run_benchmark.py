import shutil
import subprocess
import os
import time
import datetime
import sys
from threading import Thread
import re

class StartProcess(Thread):
   def __init__ (self, run_cmd, template_dir, run_dir):
      Thread.__init__(self)
      self.run_cmd = run_cmd
      self.template_dir = template_dir
      self.run_dir = run_dir
      self.status = -1
   def run(self):
      shutil.rmtree(self.run_dir, True)
      shutil.copytree(self.template_dir, self.run_dir)
      proc = subprocess.Popen(self.run_cmd, cwd=self.run_dir, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
      stdout_t, stderr_t = proc.communicate()


def compare_float(n1, n2):
    tol = 0.01
    dif_flag = True
    if n1 * n2 > 0:  #n1 and n2 have the same sign
        if abs(n1) > (1.0 + tol) * abs(n2) or abs(n1) < (1.0 - tol) * abs(n2):
            dif_flag = False
    else:
        if abs(n1 - n2) > tol:
            dif_flag = False
    return dif_flag


def comare_iobj(file1, file2):
    err_txt = ''
    try:
        with open(file1) as fin1, open(file2) as fin2:
            header1 = fin1.readline().strip()
            header2 = fin2.readline().strip()
            if not header2 or header1 != header2:
                err_txt += 'Error: Headers in .iobj files do not match\n  1: "%s"\n  2: "%s"\n' % (header1, header2)
            for line_f1 in fin1:
                line_f1 = line_f1.strip()
                line_f2 = fin2.readline().strip()
                atoms_1  = line_f1.split(',')
                atoms_2  = line_f2.split(',')
                if len(atoms_1) != len(atoms_2):
                    err_txt += 'Error: number of records do not match\n  1: "%s"\n  2: "%s"\n' % (line_f1, line_f2)
                for i1, i2 in zip(atoms_1[2:], atoms_2[2:]):
                      if not compare_float(float(i1), float(i2)):
                        err_txt += 'Error: records do not match\n  1: "%s"\n  2: "%s"\n' % (line_f1, line_f2)
    except EnvironmentError as e:
        err_txt = e.strerror
    return err_txt
    
def comare_mio(file1, file2):
    err_txt = ''
    try:
        with open(file1) as fin1, open(file2) as fin2:
            header1 = fin1.readline().strip()
            header2 = fin2.readline().strip()
            if not header2 or header1 != header2:
                err_txt += 'Error: Headers in .iobj files do not match\n  1: "%s"\n  2: "%s"\n' % (header1, header2)
            for line_f1 in fin1:
                line_f1 = line_f1.strip()
                line_f2 = fin2.readline().strip()
                if line_f1.startswith('Method of Morris for observation:') or line_f1.startswith('parameter_name, n_samples,'):
                   if not header2 or header1 != header2:
                       err_txt += 'Error: Headers in .iobj files do not match\n  1: "%s"\n  2: "%s"\n' % (header1, header2)
                else:
                    atoms_1  = line_f1.split(',')
                    atoms_2  = line_f2.split(',')
                    if len(atoms_1) != len(atoms_2):
                        err_txt += 'Error: number of records do not match\n  1: "%s"\n  2: "%s"\n' % (line_f1, line_f2)
                    for i1, i2 in zip(atoms_1[2:], atoms_1[2:]):
                        if not compare_float(float(i1), float(i2)):
                            err_txt += 'Error: records do not match\n  1: "%s"\n  2: "%s"\n' % (line_f1, line_f2)
    except EnvironmentError as e:
        err_txt = e.strerror
    return err_txt

def compare_par(file1,file2):
    err_txt = ''
    try:
        with open(file1) as fin1, open(file2) as fin2:
            header1 = fin1.readline().strip()
            header2 = fin2.readline().strip()
            if not header2 or header1 != header2:
                err_txt += 'Error: Headers in par files do not match\n  1: "%s"\n  2: "%s"\n' % (header1, header2)
            while True:
                l1,l2 = fin1.readline(),fin2.readline()
                if l1 == '' or l2 == '':
                    break
                p1 = float(l1.strip().split()[1])
                p2 = float(l2.strip().split()[1])
                if not compare_float(p1,p2):
                    err_txt += 'par {0} too different: {1}:{2}\n'.format(l1.split()[0],p1,p2)


    except EnvironmentError as e:
        err_txt += e.strerror
    return err_txt

    
def comare_sbl(file1, file2):
    err_txt = ''
    try:
        var_rex = re.compile(r'E\(.*\)\s+=\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?);\s+Var\(.*\)\s+=\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)')#\s+\(for S_i calculations\)')
        sen_rex = re.compile(r'([^,]*),\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?),\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?),\s+([-+]?[0-9]*)')
        with open(file1) as fin1, open(file2) as fin2:
            for iline1, iline2 in zip(fin1, fin2):
                iline1 = iline1.strip()
                iline2 = iline2.strip()
                m_var1 = var_rex.match(iline1)
                m_sen1 = sen_rex.match(iline1)
                if not iline1 and not iline2:
                    pass
                if iline1.startswith('Sobol Sensitivity for') or  iline2.startswith('Sobol Sensitivity for'):
                    if iline1 != iline2:
                        err_txt += 'Error: Headers in .iobj files do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                elif iline1.startswith('parameter_name,') or  iline2.startswith('parameter_name,'):
                    if iline1 != iline2:
                        err_txt += 'Error: Headers in .iobj files do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                elif m_var1:
                    m_var2 = var_rex.match(iline2)
                    if m_var2:
                        for v1, v2 in zip(m_var1.groups(), m_var2.groups()):
                             if not compare_float(float(v1), float(v2)):
                                 err_txt += 'Error: records do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                    else:
                        err_txt += 'Error: lines do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                elif m_sen1:
                    m_sen2 = sen_rex.match(iline2)
                    if m_sen2:
                        for v1, v2 in zip(m_sen1.groups()[1:3], m_sen2.groups()[1:3]):
                            if not compare_float(float(v1), float(v2)):
                                err_txt += 'Error: records do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                        if m_sen1.groups()[3] !=  m_sen2.groups()[3]:
                             err_txt += 'Error: records do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
                    else:
                        err_txt += 'Error: lines do not match\n  1: "%s"\n  2: "%s"\n' % (iline1, iline2)
    except EnvironmentError as e:
        err_txt = e.strerror
    return err_txt
    
    
def run_bm(run_dir, template_dir, exe_cmd, p_ctl, slv_exe, n_slaves):
    #assert os.path.exists(exe_cmd),"exe_cmd not found:{0}".format(exe_cmd)
    run_list = []
    #Start Master
    master_dir =  os.path.join(run_dir, 'master')
    run_cmd = r'%s %s /H :4005' % (exe_cmd, p_ctl)
    t_dir = os.path.join(run_dir, template_dir)
    print('    starting master in dir {0} with command {1}\n'.format(master_dir,run_cmd),flush=True)
    run_list.append(StartProcess(run_cmd, t_dir,  master_dir))
    run_list[-1].start()
    slave_dirs = [os.path.join(run_dir, 'slave_%d' % i) for i in range(1, n_slaves+1)]
    for idir in slave_dirs:
        run_cmd = r'%s %s /H localhost:4005' % (slv_exe, p_ctl)
        print('    starting slave in directory:{0} for command:{1}'.format(idir,run_cmd), flush=True)

        run_list.append( StartProcess(run_cmd, t_dir,  idir) )
        run_list[-1].start()
    for t in run_list:
        t.join()
    #delete slave directories
    for idir in slave_dirs:
        shutil.rmtree(idir, True)
    

if __name__ == "__main__":
    exe_cmd_pp = None
    bm_list = []
    exe_cmd_gsa = None
    bm_list_gsa = []
    f_log = open('benchmark_log.txt', 'wt')
    if sys.platform == 'win32':
       exe_cmd_pp = r'..\..\..\exe\windows\x64\Release\pestpp.exe'
       exe_cmd_gsa = r'..\..\..\exe\windows\x64\Release\gsa.exe'
       exe_cmd_opt = r'..\..\..\exe\windows\x64\Release\pestpp-opt.exe'

       #exe_cmd_pp = r'..\..\..\exe\windows\Win32\Release\pest++_32.exe'
       #exe_cmd_gsa = r'..\..\..\exe\windows\Win32\Release\gsa_32.exe'
       # run PEST++ benchmarks
       bm_list = [
        [r'.\3pg', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'.\10par_xsec', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'.\morris_1991', 'template', exe_cmd_gsa, 'pest', exe_cmd_gsa, 4, 'mio'],
        [r'.\ackley', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 10, 'iobj'],
        [r'.\box', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 10, 'iobj'],
        [r'.\kirishima', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 20, 'iobj'],

        [r'.\ishigami', 'template', exe_cmd_gsa, 'pest', exe_cmd_gsa, 10, 'sbl'],
        [r'.\opt_dewater_chance','template',exe_cmd_opt,"dewater_pest.base",exe_cmd_opt,4,'par'],
        [r'.\opt_dewater_chance', 'template', exe_cmd_opt, "dewater_pest.fosm", exe_cmd_opt, 4, 'par'],
        [r'.\opt_supply2_chance', 'template', exe_cmd_opt, "supply2_pest.base", exe_cmd_opt, 4, 'par'],
        [r'.\opt_supply2_chance', 'template', exe_cmd_opt, "supply2_pest.fosm", exe_cmd_opt, 4, 'par'],
        [r'.\opt_seawater_chance', 'template', exe_cmd_opt, "seawater_pest", exe_cmd_opt, 4, 'par']

           #[r'.\ames', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'.\tidal', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'.\hendry', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj']
        ]
    elif sys.platform == 'linux' or  sys.platform == 'linux2':
       exe_cmd_pp = r'../../../exe/linux/pestpp'
       exe_cmd_gsa = r'../../../exe/linux/gsa'
       # run PEST++ and GSA benchmarks
       bm_list = [
        [r'./stor', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'./3pg', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'./10par_xsec', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'./morris_1991', 'template_linux', exe_cmd_gsa, 'pest', exe_cmd_gsa, 4, 'mio'],
        [r'.\ackley', 'template', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'./box', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'./kirishima', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        [r'./ishigami', 'template_linux', exe_cmd_gsa, 'pest', exe_cmd_gsa, 4, 'sbl']
        #[r'./ames', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'./tidal', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj'],
        #[r'./hendry', 'template_linux', exe_cmd_pp, 'pest', exe_cmd_pp, 4, 'iobj']
        ]

    else:
        print('"%s" operating system not supported.' % sys.platform)
        print('Benchmarks will not be run')
        print()
    #print(bm_list)
    for i in bm_list:
        run_dir = i[0]
        ctl_name = i[3]
        print ('-----------------------------------------------------------------------')
        print ('starting benchmark "%s" in directory "%s"...' % (ctl_name, run_dir), flush=True)
        print(time.strftime("    %a, %d %b %Y %H:%M:%S", time.gmtime()), flush=True)
        f_log.write('-----------------------------------------------------------------------\n')
        f_log.write('starting benchmark "%s" in directory "%s"...\n' % (ctl_name, run_dir))
        f_log.write(time.strftime("    %a, %d %b %Y %H:%M:%S\n", time.gmtime()))
        start_time = time.time()
        # run benchmark
        f_log.flush()
        run_bm(*i[0:-1])
        # compare results
        err_txt = None
        if i[-1] == 'iobj':
            err_txt = comare_iobj(os.path.join(run_dir, 'baseline_opt', '%s.iobj' % ctl_name),
                              os.path.join(run_dir, 'master', '%s.iobj' % ctl_name))
        elif i[-1] == 'mio':   
            err_txt = comare_mio(os.path.join(run_dir, 'baseline_opt', '%s.mio' % ctl_name),
                              os.path.join(run_dir, 'master', '%s.mio' % ctl_name))
        elif i[-1] == 'sbl':                              
            err_txt = comare_sbl(os.path.join(run_dir, 'baseline_opt', '%s.sbl' % ctl_name),
                              os.path.join(run_dir, 'master', '%s.sbl' % ctl_name))
        elif i[-1] == "par":
            err_txt = compare_par(os.path.join(run_dir, "baseline_opt","%s.par" % ctl_name),
                                os.path.join(run_dir,"master","%s.par" % ctl_name))
        else:
            err_txt = '"%s" is an invalid comparision' % i[-1]
            
        elapsed_time = time.time() - start_time
        print(time.strftime("    %a, %d %b %Y %H:%M:%S", time.gmtime()))
        print('    elapsed time: %s' % str(datetime.timedelta(seconds=elapsed_time)))
        f_log.write(time.strftime("    %a, %d %b %Y %H:%M:%S\n", time.gmtime()))
        f_log.write('    elapsed time: %s\n' % str(datetime.timedelta(seconds=elapsed_time)))
        if  err_txt:
            print('    benchmark failed', flush=True)
            print ('  %s' % err_txt, flush=True)
            f_log.write('    benchmark failed\n')
            f_log.write ('  %s\n' % err_txt)
        else:
            print('    benchmark passed', flush=True)
            f_log.write('    benchmark passed\n')
            shutil.rmtree(os.path.join(run_dir, 'master'), True)
        print(flush=True)
        f_log.write('\n')
        f_log.flush()
    f_log.close()

        
 
