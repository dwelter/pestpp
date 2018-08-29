import numpy as np
from numpy.ctypeslib import ndpointer 
import ctypes
import os
import copy

class RunpyError(Exception):
    def __init__(self, message):
        self.message = message



class runpy:
    def __init__(self):
        self.runlib = ctypes.cdll.LoadLibrary(r'wrappers.dll')
        self.obj = ctypes.c_void_p()
        print('runpy started')

    def string2c(self,str):
       b_str = str.encode('utf-8')
       return b_str
   
    def list_strings2c(self,L):
       arr = (ctypes.c_char_p * len(L))()
       for i in range(0, len(arr)):
          arr[i] = self.string2c(L[i])
       return arr

    def create_serial(self, comline, tpl, inp, ins, out, stor_file, rundir, n_fail_max):
        func = self.runlib.rmic_create_serial
        func.restype = ctypes.c_void_p
        comline_c = self.list_strings2c(comline)
        tpl_c = self.list_strings2c(tpl)
        inp_c = self.list_strings2c(inp)
        ins_c = self.list_strings2c(ins)
        out_c = self.list_strings2c(out)
        stor_file_c = self.string2c(stor_file)
        rundir_c = self.string2c(rundir)
        self.obj.value = func(comline_c, len(comline),
              tpl_c, len(tpl), 
              inp_c, len(inp),
              ins_c, len(ins),
              out_c, len(out), 
              stor_file_c, rundir_c, n_fail_max)
        
    def create_panther(self, stor_file, port, info_file, 
                    n_fail_max, overdue_resched_fac, overdue_giveup_fac):
        func = self.runlib.rmic_create_panther
        func.restype = ctypes.c_void_p
        stor_file_c = self.string2c(stor_file)
        port_c = self.string2c(port)
        info_file_c = self.string2c(info_file)
        self.obj.value = func(stor_file_c, port_c, info_file_c, n_fail_max, 
                  ctypes.c_double(overdue_resched_fac), ctypes.c_double(overdue_giveup_fac))


    def initialize(self, par_names, obs_names):
        self.par_names = copy.deepcopy(par_names)
        self.obs_names = copy.deepcopy(obs_names)
        func = self.runlib.rmic_initialize
        func.restype = ctypes.c_int
        par_names_c = self.list_strings2c(self.par_names)
        obs_names_c = self.list_strings2c(self.obs_names)
        err = func(self.obj, par_names_c, len(self.par_names), obs_names_c, len(self.obs_names))
        if err != 0:
            raise RunpyError('runpy error in function runpy.intialize(): error number = %d' % err)

    def reinitialize(self):
        func = self.runlib.rmic_initialize
        func.restype = ctypes.c_int
        err = func(self.obj)
        if err != 0:
            raise RunpyError('runpy error in function runpy.reinitialize(): error number = %d' % err)

    def add_run(self, par_data, model_exe_index):
        func = self.runlib.rmic_add_run
        func.restype = ctypes.c_int
        func.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int,ctypes.POINTER(ctypes.c_int)]
        n_par = ctypes.c_int(par_data.size)
        run_id =  ctypes.c_int(-1)
        err = func(self.obj, par_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                         n_par, model_exe_index, ctypes.byref(run_id))
        if err != 0:
            raise RunpyError('runpy error in function runpy.add_run(): error number = %d' % err)
        return run_id.value

    def run(self):
        func = self.runlib.rmic_run
        func.restype = ctypes.c_int
        func.argtypes = [ctypes.c_void_p]
        err = func(self.obj)
        if err != 0:
            raise RunpyError('runpy error in function runpy.run(): error number = %d' % err)
        return err

    def get_run(self, run_id):
        func = self.runlib.rmic_get_run
        func.restype = ctypes.c_int
        func.argtypes = [ctypes.c_void_p, ctypes.c_int, 
                 ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                 ctypes.POINTER(ctypes.c_double), ctypes.c_int]
        run_id_c = ctypes.c_int(run_id)
        n_par_c = ctypes.c_int(len(self.par_names))
        n_obs_c = ctypes.c_int(len(self.obs_names))
        par_data = np.zeros(len(self.par_names), dtype=np.double)
        obs_data = np.zeros(len(self.obs_names), dtype=np.double)
        err = func(self.obj, run_id_c, 
            par_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_par_c,
            obs_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_obs_c)
        if err != 0:
            #raise RunpyError('runpy error in function runpy.get_run() with id = %d: error number = %d' % (run_id, err))
            #par_data = None
            obs_data = None
        return par_data, obs_data

    def get_num_failed(self):
        func = self.runlib.rmic_get_num_failed_runs
        func.restype = ctypes.c_int
        nfail = ctypes.c_int(-999)
        err = func(self.obj, ctypes.byref(nfail))
        if err != 0:
            raise RunpyError('runpy error in function runpy.get_num_failed() with id = %d: error number = %d' % (run_id, err))
        return nfail.value

    def get_failed_run_ids(self):
        nfail = self.get_num_failed()
        fail_list = []
        func = self.runlib.rmic_get_failed_runs_n
        func.restype = ctypes.c_int
        func.argtypes = [ ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int ]
        fail_array = np.ones((nfail,), dtype=np.int)
        err = func(self.obj, fail_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), ctypes.c_int(nfail))
        fail_list = fail_array.tolist()
        if err != 0:
            raise RunpyError('runpy error in function runpy.get_failed_run_ids() with id = %d: error number = %d' % (run_id, err))
        return fail_list

    def __del__(self):
        func = self.runlib.rmic_delete
        func.restype = ctypes.c_int
        func.argtypes = [ctypes.c_void_p]
        err = func(self.obj)
        if err != 0:
            raise RunpyError('runpy error in function runpy.__del__(): error number = %d' % err)



def test_serial():
    comline = [r'.\exe\storage1.exe']
    tpl = [r'.\misc\input.tpl']
    inp = [r'input.dat']
    ins = [r'.\misc\output.ins']
    out = [r'output.dat']
    stor_file = 'tmp_run_data.bin'
    #rundir = r'..\..\benchmarks\stor\test'
    rundir = r'C:\Users\Dave\Documents\GitHub\pestpp\benchmarks\stor\test'

    par_names = ['recharge', 'cond', 'scoeff']
    obs_names = ['head1', 'head2', 'head3', 'head4',
                 'head5', 'head6', 'head7', 'head8',
                 'head9', 'head10', 'head11', 'head12',
                 'head13', 'head14', 'head15', 'head16']

    run_mngr = runpy()
    run_mngr.create_serial(comline, tpl, inp, ins, out, stor_file, rundir, 1)
    run_mngr.initialize(par_names, obs_names)
    #add 3 runs
    nruns = 0
    for i in range(0, 3):
        data = np.array([0.1 + .02*i, 0.005, 0.05], dtype=np.double)
        run_mngr.add_run(data, 1)
        nruns += 1
    data = np.array([-9.0e99, 0.005, 0.05], dtype=np.double)
    run_mngr.add_run(data, 1)
    nruns += 1

    err = run_mngr.run()

    failed_runs = run_mngr.get_failed_run_ids()
    print('failed runs:', failed_runs)
    for i in range(0, nruns):
        print('Results for run: %d' % i)
        par, obs = run_mngr.get_run(i)
        print(par)
        print(obs)
        print('')


def test_panther():
    stor_file = 'tmp_run_data.bin'
    port = '4005'
    out_file = 'panther_test.out'
    par_names = ['recharge', 'cond', 'scoeff']
    obs_names = ['head1', 'head2', 'head3', 'head4',
                 'head5', 'head6', 'head7', 'head8',
                 'head9', 'head10', 'head11', 'head12',
                 'head13', 'head14', 'head15', 'head16']

    run_mngr = runpy()
    run_mngr.create_panther(stor_file, port, out_file, 3, 1.15, 100.0)
    run_mngr.initialize(par_names, obs_names)
    #add 3 runs
    nruns = 3
    for i in range(0, nruns):
        data = np.array([0.1 * .02*nruns, 0.005, 0.05], dtype=np.double)
        run_mngr.add_run(data,1)

    run_mngr.run()

    for i in range(0, nruns):
        print('Results for run: %d' % i)
        par, obs = run_mngr.get_run(i)
        print(par)
        print(obs)
        print('')
    pass

if __name__=="__main__":
    print('Starting serial test')
    test_serial() 
    print('Starting PANTHER test')
    test_panther()