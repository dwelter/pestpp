import os
import shutil
import numpy as np
import pandas as pd
import flopy
import pyemu

pyemu_dir = "es_pyemu"
pestpp_dir = "master"

def prep():
    if os.path.exists(pyemu_dir):
        shutil.rmtree(pyemu_dir)
    if os.path.exists(pestpp_dir):
        shutil.rmtree(pestpp_dir)

    shutil.copytree("template",pyemu_dir)
    shutil.copytree("template", pestpp_dir)
    num_reals = 15
    pst = pyemu.Pst(os.path.join(pyemu_dir,"pest.pst"))
    pst.control_data.noptmax = 20
    pst.pestpp_options["ies_parameter_csv"] = "par.csv"
    pst.pestpp_options["ies_observation_csv"] = "obs.csv"
    pst.pestpp_options["ies_obs_restart_csv"] = "restart_obs.csv"
    pst.pestpp_options["ies_use_approx"] = "true"
    pst.pestpp_options["ies_use_prior_scaling"] = "false"
    pst.pestpp_options["ies_num_reals"] = "{0}".format(num_reals)
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.svd_data.eigthresh = 1.0e-5
    pst.svd_data.maxsing = 20
    #pst.observation_data.loc[pst.nnz_obs_names,"weight"] /= 10.0
    pst.write(os.path.join(pyemu_dir,"pest.pst"))
    pst.write(os.path.join(pestpp_dir, "pest.pst"))
    #dia_parcov = pyemu.Cov.from_parameter_data(pst,sigma_range=6.0)
    parcov = pyemu.Cov.from_parameter_data(pst)
    print(parcov.as_2d)
    parcov = pyemu.Cov(parcov.as_2d,names=parcov.row_names)
    print(parcov.x)
    parcov.to_ascii(os.path.join(pestpp_dir,"prior.cov"))

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,parcov,num_reals=num_reals,use_homegrown=True)
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst,num_reals=num_reals)
    oe.to_csv(os.path.join(pyemu_dir,"obs.csv"))
    oe.to_csv(os.path.join(pestpp_dir, "obs.csv"))
    oe.to_binary(os.path.join(pestpp_dir, "obs.jcb"))

    pe.to_csv(os.path.join(pyemu_dir,"sweep_in.csv"))
    pe.to_csv(os.path.join(pyemu_dir,"par.csv"))
    pe.to_csv(os.path.join(pestpp_dir, "par.csv"))
    pe.to_binary(os.path.join(pestpp_dir, "par.jcb"))

    pyemu.helpers.start_slaves(pyemu_dir,"sweep","pest.pst", num_slaves=5,master_dir="sweep_master",cleanup=False)
    df = pd.read_csv(os.path.join("sweep_master","sweep_out.csv"))
    df.to_csv(os.path.join(pyemu_dir,"restart_obs.csv"))
    df = df.loc[:,pst.observation_data.obsnme.apply(str.upper)]
    df.to_csv(os.path.join(pestpp_dir,"restart_obs.csv"))
    oe = pyemu.ObservationEnsemble.from_dataframe(df=df,pst=pst)
    oe.to_binary(os.path.join(pestpp_dir, "restart_obs.jcb"))

def ies():
    os.chdir(pyemu_dir)
    es = pyemu.EnsembleSmoother("pest.pst",verbose="ies.log",save_mats=True,
                                slave_dir=os.path.join("..","template"),num_slaves=5)
    es.initialize(parensemble="par.csv",obsensemble="obs.csv",restart_obsensemble="restart_obs.csv")
    for i in range(es.pst.control_data.noptmax):
        es.update(use_approx=False)
    #es.update(lambda_mults=[0.1,1.0,# 10.0],run_subset=10)

    #es.update(lambda_mults=[0.1,1.0,10.0],run_subset=10)
    os.chdir('..')

def test():
    pst = pyemu.Pst(os.path.join(pestpp_dir,"pest.pst"))
    m = pyemu.ObservationEnsemble.from_binary(pst,os.path.join(pestpp_dir, "restart_obs.jcb"))
    print(m)

def check():
    pst = pyemu.Pst(os.path.join(pyemu_dir, "pest.pst"))
    parcov = pyemu.Cov.from_parameter_data(pst)
    df = pd.read_csv(os.path.join("master","pest.0.par.csv"))
    df = df.apply(np.log10)
    print(df.std())
    print(df.mean())

def draw_invest():
    pst = pyemu.Pst(os.path.join(pyemu_dir, "pest.pst"))
    pst.control_data.noptmax = 20
    try:
        pst.pestpp_options.pop("ies_parameter_csv")
        pst.pestpp_options.pop("ies_observation_csv")
        pst.pestpp_options.pop("ies_obs_restart_csv")
    except:
        pass
    pst.pestpp_options["ies_use_approx"] = "true"
    pst.pestpp_options["ies_use_prior_scaling"] = "false"
    pst.pestpp_options["ies_num_reals"] = "10";
    pst.svd_data.eigthresh = 1.0e-5
    pst.svd_data.maxsing = 20
    # pst.observation_data.loc[pst.nnz_obs_names,"weight"] /= 10.0
    pst.write(os.path.join(pyemu_dir, "pest.pst"))
    pst.write(os.path.join(pestpp_dir, "pest.pst"))
    # dia_parcov = pyemu.Cov.from_parameter_data(pst,sigma_range=6.0)
    parcov = pyemu.Cov.from_parameter_data(pst)
    num_reals = 10
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, parcov, num_reals=num_reals, use_homegrown=True)
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst, num_reals=num_reals)
if __name__ == "__main__":
    #prep()
    #ies()
    #test()
    pyemu.helpers.start_slaves("template","sweep.exe","pest.pst",num_slaves=5)
    #draw_invest()
    #check()