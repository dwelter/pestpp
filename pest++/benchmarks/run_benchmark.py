import shutil
import subprocess
import os
import time
from threading import Thread


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
    
def run_bm(run_dir, template_dir, exe_cmd, p_ctl, n_slaves):
    run_list = []
    #Start Master
    master_dir =  r'%s\master' % run_dir
    run_cmd = r'%s %s /H :4005' % (exe_cmd, p_ctl)
    t_dir = '%s\%s' % (run_dir, template_dir)
    run_list.append(StartProcess(run_cmd, t_dir,  master_dir))
    run_list[-1].start()
    slave_dirs = [r"%s\slave_%d\\" % (run_dir, i) for i in range(1, n_slaves+1)]
    for idir in slave_dirs:
        run_cmd = r'%s %s /H localhost:4005' % (exe_cmd, p_ctl)
        run_list.append( StartProcess(run_cmd, t_dir,  idir) )
        run_list[-1].start()
    for t in run_list:
        t.join()
    #delete slave directories
    for idir in slave_dirs:
        shutil.rmtree(idir, True)
    

if __name__ == "__main__":
    exe_cmd = r'..\..\..\x64\Release\pest++.exe'
    # run PEST++ benchmarks
    bm_list = [
        [r'.\stor', 'template', exe_cmd, 'pest', 4],
        [r'.\3pg', 'template', exe_cmd, 'pest', 4],
        [r'.\10par_xsec', 'template', exe_cmd, 'pest', 4],
        [r'.\box', 'template', exe_cmd, 'pest', 4],
        [r'.\kirishima', 'template', exe_cmd, 'pest', 4],
        [r'.\ishigami', 'template', exe_cmd, 'pest', 4]
        #[r'.\ames', 'template', exe_cmd, 'pest', 4],
        #[r'.\tidal', 'template', exe_cmd, 'pest', 4],
        #[r'.\hendry', 'template', exe_cmd, 'pest', 4]
    ]
    for i in bm_list:
        run_dir = i[0]
        ctl_name = i[3]
        print ('starting benchmark in directory "%s"...' % run_dir, flush=True)
        # run benchmark
        run_bm(*i)
        # compare results
        err_txt = comare_iobj(r'%s\baseline_opt\%s.iobj' % (run_dir, ctl_name), r'%s\master\%s.iobj' % ( run_dir, ctl_name))
        if  err_txt:
            print('  benchmark failed', flush=True)
            print ('  %s' % err_txt, flush=True)
        else:
            print('  benchmark passed', flush=True)
        print(flush=True)
        #clean up
        #delete the master if the benchmarch passed
        if not err_txt:
            shutil.rmtree(r'%s\master' % run_dir , True)
        
        
    # run GSA++ benchmarks
    exe_cmd = r'..\..\..\x64\Release\gsa.exe'
    bm_list = [
            [r'.\morris_1991', 'template', exe_cmd, 'pest', 4]
    ]
    for i in bm_list:
        run_dir = i[0]
        ctl_name = i[3]
        print ('starting benchmark in directory "%s"...' % run_dir, flush=True)
        # run benchmark
        run_bm(*i)
        # compare results
        err_txt = comare_mio(r'%s\baseline_opt\%s.mio' % (run_dir, ctl_name), r'%s\master\%s.mio' % ( run_dir, ctl_name))
        if  err_txt:
            print('  benchmark failed', flush=True)
            print ('  %s' % err_txt, flush=True)
        else:
            print('  benchmark passed', flush=True)
        print(flush=True)
        #clean up
        #delete the master if the benchmarch passed
        if not err_txt:
            shutil.rmtree(r'%s\master' % run_dir , True)