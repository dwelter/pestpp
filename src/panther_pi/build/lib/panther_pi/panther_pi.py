import numpy as np
from numpy.ctypeslib import ndpointer 
import ctypes as ct
import os
import copy

class PantherError(Exception):
    def __init__(self, message):
        self.message = message

class panther:
    def __init__(self):
        dll_file = os.path.join( os.path.dirname(__file__), 'wrappers.dll')
        self.runlib = ct.cdll.LoadLibrary(dll_file)
        self.obj = ct.c_void_p()

    def string2c(self,str):
       b_str = str.encode('utf-8')
       return b_str
   
    def list_strings2c(self,L):
       arr = (ct.c_char_p * len(L))()
       for i in range(0, len(arr)):
          arr[i] = self.string2c(L[i])
       return arr
   
    def err_msg(self):
        func = self.runlib.rmic_err_msg
        func.argtypes = []
        func.restype = ct.c_char_p
        c_err_msg = func()
        return c_err_msg.decode('utf-8')

    def create_serial(self, comline, tpl, inp, ins, out, stor_file, rundir, n_fail_max):
        func = self.runlib.rmic_create_serial
        func.argtypes = [ct.POINTER(ct.c_char_p), ct.c_int,
                         ct.POINTER(ct.c_char_p), ct.c_int,
                         ct.POINTER(ct.c_char_p), ct.c_int,
                         ct.POINTER(ct.c_char_p), ct.c_int,
                         ct.POINTER(ct.c_char_p), ct.c_int,
                         ct.c_char_p, ct.c_char_p, ct.c_int]
        func.restype = ct.c_void_p
        comline_c = self.list_strings2c(comline)
        tpl_c = self.list_strings2c(tpl)
        inp_c = self.list_strings2c(inp)
        ins_c = self.list_strings2c(ins)
        out_c = self.list_strings2c(out)
        stor_file_c = self.string2c(stor_file)
        rundir_c = self.string2c(rundir)
        self.obj.value = func(comline_c, ct.c_int(len(comline)),
              tpl_c, ct.c_int(len(tpl)), 
              inp_c, ct.c_int(len(inp)),
              ins_c, ct.c_int(len(ins)),
              out_c, ct.c_int(len(out)), 
              stor_file_c, rundir_c, ct.c_int(n_fail_max))
        
    def create_panther(self, stor_file, port, info_file, 
                    n_fail_max, overdue_resched_fac, overdue_giveup_fac):
        func = self.runlib.rmic_create_panther
        func.argtypes = [ct.c_char_p, ct.c_char_p, ct.c_char_p,
                         ct.c_int, ct.c_double, ct.c_double]
        func.restype = ct.c_void_p
        stor_file_c = self.string2c(stor_file)
        port_c = self.string2c(port)
        info_file_c = self.string2c(info_file)
        self.obj.value = func(stor_file_c, port_c, info_file_c, ct.c_int(n_fail_max), 
                  ct.c_double(overdue_resched_fac), ct.c_double(overdue_giveup_fac))

    def initialize(self, par_names, obs_names):
        self.par_names = copy.deepcopy(par_names)
        self.obs_names = copy.deepcopy(obs_names)
        func = self.runlib.rmic_initialize
        func.argtypes = [ ct.c_void_p,
                          ct.POINTER(ct.c_char_p), ct.c_int,
                          ct.POINTER(ct.c_char_p), ct.c_int ]
        func.restype = ct.c_int
        par_names_c = self.list_strings2c(self.par_names)
        obs_names_c = self.list_strings2c(self.obs_names)
        err = func(self.obj, par_names_c, ct.c_int(len(self.par_names)), obs_names_c, ct.c_int(len(self.obs_names)))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.intialize(): error number = %d\n %s' % (err, self.err_msg()))

    def reinitialize(self):
        func = self.runlib.rmic_initialize
        func.argtypes = [ct.c_void_p]
        func.restype = ct.c_int
        err = func(self.obj)
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.reinitialize(): error number = %d\n %s' % (err, self.err_msg()))
            
    def initialize_restart(self, store_filename):
        func = self.runlib.intiialize_restart
        func.argtypes = [ct.c_void_p, ct.c_char_p]
        func.restype = ct.c_int
        stor_file_c = self.string2c(stor_file)
        err = func(self.obj, stor_file_c)
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.intialize(): error number = %d\n %s' % (err, self.err_msg()))

    def add_run(self, par_data, model_exe_index):
        func = self.runlib.rmic_add_run
        func.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), ct.c_int, ct.c_int,ct.POINTER(ct.c_int)]
        func.restype = ct.c_int
        n_par = ct.c_int(par_data.size)
        run_id =  ct.c_int(-1)
        err = func(self.obj, par_data.ctypes.data_as(ct.POINTER(ct.c_double)),
                         n_par, model_exe_index, ct.byref(run_id))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.add_run(): error number = %d\n %s' % (err, self.err_msg()))
        return run_id.value
    
    def add_run_with_info(self, par_data, model_exe_index, info_txt, info_value):
        func = self.runlib.rmic_add_run_with_info
        func.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double), ct.c_int, ct.c_int,
                         ct.c_char_p, ct.c_double, ct.POINTER(ct.c_int)]
        func.restype = ct.c_int
        n_par = ct.c_int(par_data.size)
        model_exe_index_c = ct.c_int(model_exe_index)
        info_txt_c = self.string2c(info_txt)
        info_value_c = ct.c_double(info_value)
        run_id =  ct.c_int(-1)
        err = func(self.obj, par_data.ctypes.data_as(ct.POINTER(ct.c_double)),
                         n_par, model_exe_index_c, info_txt_c, info_value_c, ct.byref(run_id))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.add_run(): error number = %d\n %s' % (err, self.err_msg()))
        return run_id.value

    def run(self):
        func = self.runlib.rmic_run
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p]
        err = func(self.obj)
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.run(): error number = %d\n %s' % (err, self.err_msg()))
    
    def run_until(self, condition, no_ops, time_sec):
        return_cond = 0
        func = self.runlib.rmic_run
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, ct.POINTER(ct.c_int)]
        return_cond_c = ct.c_int(0)
        err = func(self.obj, ct.c_int(no_ops), ct.c_double(time_sec), ct.byref(return_cond_c))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.run_until(): error number = %d\n %s' % (err, self.err_msg()))
        return return_cond
    
    def cancel_run(self, run_id):
        func = self.runlib.rmic_cancel_run
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int]
        err = func(self.obj, ct.c_int(run_id))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.cancel_run(): error number = %d\n %s' % (err, self.err_msg()))

    def get_run_status(self, run_id):
        func = self.runlib.rmic_get_run_status
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)]
        run_status =  ct.c_int(-1)
        max_runtime = ct.c_double(-1)
        n_concurrent_runs = ct.c_int(-1)
        err = func(self.obj, ct.c_int(run_id), ct.byref(run_status), ct.byref(max_runtime), ct.byref(n_concurrent_runs))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.get_run_status(): error number = %d\n %s' % (err, self.err_msg()))
        return run_status.value, max_runtime.value, n_concurrent_runs.value
        
    def get_run_info(self, run_id):
        func = self.runlib.rmic_get_run_info
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_char_p , ct.POINTER(ct.c_int), ct.POINTER(ct.c_double)]
        run_status =  ct.c_int(-1)
        model_exe_index = ct.c_int(-1)
        info_txt_c = self.string2c(''*41)
        info_txt_c_len = ct.c_int(41)
        n_concurrent_runs = ct.c_int(-1)
        info_value = ct.c_double(-1)
        err = func(self.obj, ct.c_int(run_id), ct.byref(run_status), ct.byref(model_exe_index), info_txt_c, info_txt_c_len, ct.byref(info_value))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.get_run_status(): error number = %d\n %s' % (err, self.err_msg()))
        return run_status.value, max_runtime.value, n_concurrent_runs.value   

    def get_run(self, run_id):
        func = self.runlib.rmic_get_run
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int, 
                 ct.POINTER(ct.c_double), ct.c_int,
                 ct.POINTER(ct.c_double), ct.c_int]
        run_id_c = ct.c_int(run_id)
        n_par_c = ct.c_int(len(self.par_names))
        n_obs_c = ct.c_int(len(self.obs_names))
        par_data = np.zeros(len(self.par_names), dtype=np.double)
        obs_data = np.zeros(len(self.obs_names), dtype=np.double)
        err = func(self.obj, run_id_c, 
            par_data.ctypes.data_as(ct.POINTER(ct.c_double)), n_par_c,
            obs_data.ctypes.data_as(ct.POINTER(ct.c_double)), n_obs_c)
        #err = err_c.value
        if err != 0:
             run_ok = False
             obs_data = None
             #raise PantherError('panther error in function panther.get_run() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
        return par_data, obs_data
    
    def get_run_with_info(run_id):
        func = self.runlib.rmic_get_run_with_info
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p, ct.c_int, 
                 ct.POINTER(ct.c_double), ct.c_int,
                 ct.POINTER(ct.c_double), ct.c_int]
        run_id_c = ct.c_int(run_id)
        n_par_c = ct.c_int(len(self.par_names))
        n_obs_c = ct.c_int(len(self.obs_names))
        par_data = np.zeros(len(self.par_names), dtype=np.double)
        obs_data = np.zeros(len(self.obs_names), dtype=np.double)
        info_txt_c = self.string2c(info_txt)
        info_value_c = ct.c_double(info_value)
        
        err = func(self.obj, run_id_c, 
            par_data.ctypes.data_as(ct.POINTER(ct.c_double)), n_par_c,
            obs_data.ctypes.data_as(ct.POINTER(ct.c_double)), n_obs_c,
            ct.byref(info_txt_c), ct,byref(info_value_c))
        
        info_txt = info_txt_c.value
        info_value = info_value_c.value
        #err = err_c.value
        run_ok = True
        if err != 0:
             run_ok = False
             #raise PantherError('panther error in function panther.get_run_with_info() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
             obs_data = None
        return par_data, obs_data, info_txt, info_value

    def get_n_failed_runs(self):
        func = self.runlib.rmic_get_n_failed_runs
        func.restype = ct.c_int
        nfail = ct.c_int(-999)
        err = func(self.obj, ct.byref(nfail))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.get_num_failed() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
        return nfail.value

    def get_failed_run_ids(self):
        nfail = self.get_n_failed_runs()
        fail_list = []
        func = self.runlib.rmic_get_failed_run_ids
        func.restype = ct.c_int
        func.argtypes = [ ct.c_void_p, ct.POINTER(ct.c_int), ct.c_int ]
        fail_array = np.ones((nfail,), dtype=np.int)
        err = func(self.obj, fail_array.ctypes.data_as(ct.POINTER(ct.c_int)), ct.c_int(nfail))
        #err = err_c.value
        fail_list = fail_array.tolist()
        if err != 0:
            raise PantherError('panther error in function panther.get_failed_run_ids() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
        return fail_list

    def get_n_cur_runs():
        func = self.runlib.rmic_get_n_cur_runs
        func.restype = ct.c_int
        ncur = ct.c_int(-999)
        err = func(self.obj, ct.byref(ncur))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.get_n_cur_runs() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
        return ncur.value

    
    def get_n_total_runs():
        func = self.runlib.rmic_get_n_total_runs
        func.restype = ct.c_int
        nruns = ct.c_int(-999)
        err = func(self.obj, ct.byref(nruns))
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.get_n_total_runs() with id = %d: error number = %d\n%s' % (run_id, err, self.err_msg()))
        return nruns.value
    

    def __del__(self):
        func = self.runlib.rmic_delete
        func.restype = ct.c_int
        func.argtypes = [ct.c_void_p]
        err = func(self.obj)
        #err = err_c.value
        if err != 0:
            raise PantherError('panther error in function panther.__del__(): error number = %d\n %s' % (err, self.err_msg()))



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

    run_mngr = panther()
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

    run_mngr = panther()
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

if __name__=="__main__":
    print('Starting serial test')
    test_serial() 
    print('Starting PANTHER test')
    test_panther()