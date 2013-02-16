import numpy as np
from numpy.ctypeslib import ndpointer 
import ctypes
import os
print os.getcwd()
import copy


#convert a list of python strings to an array of ctypes char_p
def call_c(L):
    arr = (ctypes.c_char_p * len(L))()
    arr[:] = L
    return arr

class runpy:
    def __init__(self):
        self.runlib = ctypes.cdll.LoadLibrary('.\\yamr')
        self.obj = ctypes.c_void_p()

    def call_c(self,L):
        arr = (ctypes.c_char_p * len(L))()
        arr[:] = L
        return arr

    def create_serial(self, comline, tpl, inp, ins, out, stor_file, rundir):
        func = self.runlib.rmic_create_serial
        func.restype = ctypes.c_void_p
        comline_c = self.call_c(comline)
        tpl_c = self.call_c(tpl)
        inp_c = self.call_c(inp)
        ins_c = self.call_c(ins)
        out_c = self.call_c(out)
        self.obj = func(comline_c, len(comline),
              tpl_c, len(tpl), 
              inp_c, len(inp),
              ins_c, len(ins),
              out_c, len(out), 
              ctypes.c_char_p(stor_file), ctypes.c_char_p(rundir))


    def initialize(self, par_names, obs_names):
        self.par_names = copy.deepcopy(par_names)
        self.obs_names = copy.deepcopy(obs_names)
        func = self.runlib.rmic_initialize
        par_names_c = self.call_c(self.par_names)
        obs_names_c = self.call_c(self.obs_names)
        func(self.obj, par_names_c, len(self.par_names), obs_names_c, len(self.obs_names))

    def add_run(self, par_data):
        func = self.runlib.rmic_add_run
        func.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
        n_par = ctypes.c_int(par_data.size)
        run_id =  ctypes.c_int(-1)
        func(self.obj, par_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                         n_par, ctypes.byref(run_id))
        return run_id.value

    def run(self):
        func = self.runlib.rmic_run
        func(self.obj)

    def get_run(self, run_id):
        func = self.runlib.rmic_get_run
        func.argtypes = [ctypes.c_void_p, ctypes.c_int, 
                 ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                 ctypes.POINTER(ctypes.c_double), ctypes.c_int]
        run_id_c = ctypes.c_int(run_id)
        n_par_c = ctypes.c_int(len(self.par_names))
        n_obs_c = ctypes.c_int(len(self.obs_names))
        par_data = np.zeros(len(self.par_names), dtype=np.double)
        obs_data = np.zeros(len(self.obs_names), dtype=np.double)
        func(self.obj, run_id_c, 
            par_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_par_c,
            obs_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_obs_c)
        return par_data, obs_data



def test():
    comline = ['storage1.exe']
    tpl = ['input.tpl']
    inp = ['input.dat']
    ins = ['output.ins']
    out = ['output.dat']
    stor_file = 'tmp_run_data.bin'
    rundir = r'C:\Users\Dave\Documents\Visual Studio 2012\Projects\pestpp\pest++\benchmarks\stor'

    par_names = ['recharge', 'cond', 'scoeff']
    obs_names = ['head1', 'head2', 'head3', 'head4',
                 'head5', 'head6', 'head7', 'head8',
                 'head9', 'head10', 'head11', 'head12',
                 'head13', 'head14', 'head15', 'head16']

    data = np.array([0.1, 0.005, 0.05], dtype=np.double)

    run_mngr = runpy()
    run_mngr.create_serial(comline, tpl, inp, ins, out, stor_file, rundir)
    run_mngr.initialize(par_names, obs_names)
    run_mngr.add_run(data)
    run_mngr.run()
    par, obs = run_mngr.get_run(0)


test() 