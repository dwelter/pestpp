import os
import shutil
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

    #os.chdir(os.path.join("smoother","freyberg"))

    #if not os.path.exists("freyberg.xy"):

    ml = flopy.modflow.Modflow.load("freyberg.nam",model_ws=pyemu_dir,
                                    load_only=[])
    #xy = pd.DataFrame([(x,y) for x,y in zip(ml.sr.xcentergrid.flatten(),ml.sr.ycentergrid.flatten())],
    #                  columns=['x','y'])
    x,y,names = [],[],[]
    for i in range(ml.nrow):
        for j in range(ml.ncol ):
            if ml.bas6.ibound.array[0,i,j] <1:
                continue
            names.append("hkr{0:02d}c{1:02d}".format(i,j))
            x.append(ml.sr.xcentergrid[i,j])
            y.append(ml.sr.ycentergrid[i,j])
        #xy.loc[:,"name"] = names
        #xy.to_csv("freyberg.xy")
    #else:
        #xy = pd.read_csv("freyberg.xy")

    pst = pyemu.Pst(os.path.join(pyemu_dir,"pest.pst"))
    pst.control_data.noptmax = 6
    pst.pestpp_options = {}
    #pst.pestpp_options["ies_parameter_csv"] = "par.csv"
    #pst.pestpp_options["ies_observation_csv"] = "obs.csv"
    #pst.pestpp_options["ies_obs_restart_csv"] = "restart_obs.csv"
    pst.pestpp_options["parcov_filename"] = "freyberg_prior.jcb"
    pst.pestpp_options["ies_use_approx"] = "true"
    pst.svd_data.eigthresh = 1.0e-4
    pst.svd_data.maxsing = 1000000
    #pst.observation_data.loc[pst.nnz_obs_names,"weight"] /= 10.0
    pst.write(os.path.join(pyemu_dir,"pest.pst"))
    pst.write(os.path.join(pestpp_dir, "pest.pst"))
    #dia_parcov = pyemu.Cov.from_parameter_data(pst,sigma_range=6.0)
    par = pst.parameter_data
    hk_names = par.loc[par.parnme.apply(lambda x: x.startswith("hk")),'parnme']
    hk_par = par.loc[hk_names,:]
    hk_par.loc[:,"i"] = hk_par.parnme.apply(lambda x: int(x[3:5]))
    hk_par.loc[:, "j"] = hk_par.parnme.apply(lambda x: int(x[-2:]))
    hk_par.loc[:,"x"] = hk_par.apply(lambda x: ml.sr.xcentergrid[x.i,x.j],axis=1)
    hk_par.loc[:, "y"] = hk_par.apply(lambda x: ml.sr.ycentergrid[x.i, x.j], axis=1)


    gs = pyemu.utils.geostats.read_struct_file(os.path.join(pyemu_dir,"structure.dat"))
    parcov = pyemu.helpers.geostatistical_prior_builder(pst,struct_dict={gs:hk_par},sigma_range=6)

    parcov.to_binary(os.path.join(pyemu_dir,"freyberg_prior.jcb"))
    parcov.to_binary(os.path.join(pestpp_dir,"freyberg_prior.jcb"))

    obscov = pyemu.Cov.from_observation_data(pst)
    num_reals = 15
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,parcov,num_reals=num_reals,use_homegrown=True)
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst,num_reals=num_reals)
    oe.to_csv(os.path.join(pyemu_dir,"obs.csv"))
    oe.to_csv(os.path.join(pestpp_dir, "obs.csv"))

    pe.to_csv(os.path.join(pyemu_dir,"sweep_in.csv"))
    pe.to_csv(os.path.join(pyemu_dir,"par.csv"))
    pe.to_csv(os.path.join(pestpp_dir, "par.csv"))

    pyemu.helpers.start_slaves(pyemu_dir,"sweep","pest.pst", num_slaves=10,master_dir="sweep_master",cleanup=False)
    df = pd.read_csv(os.path.join("sweep_master","sweep_out.csv"))
    df.to_csv(os.path.join(pyemu_dir,"restart_obs.csv"))
    df = df.loc[:,pst.observation_data.obsnme.apply(str.upper)]
    df.to_csv(os.path.join(pestpp_dir,"restart_obs.csv"))
def ies():
    os.chdir(pyemu_dir)
    parcov = pyemu.Cov.from_binary("freyberg_prior.jcb")
    es = pyemu.EnsembleSmoother("pest.pst",parcov=parcov,verbose="ies.log",save_mats=True,
                                slave_dir=os.path.join("..","template"),num_slaves=5)
    es.initialize(parensemble="par.csv",obsensemble="obs.csv",restart_obsensemble="restart_obs.csv")
    for i in range(es.pst.control_data.noptmax):
        es.update(use_approx=True)
    #es.update(lambda_mults=[0.1,1.0,# 10.0],run_subset=10)

    #es.update(lambda_mults=[0.1,1.0,10.0],run_subset=10)
    os.chdir('..')
if __name__ == "__main__":
    #prep()
    #ies()
    pyemu.helpers.start_slaves("template","sweep.exe","pest.pst",num_slaves=10)