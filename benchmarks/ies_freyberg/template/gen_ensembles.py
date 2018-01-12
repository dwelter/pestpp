import os
import flopy
import pyemu

def prep():
    import os
    import pandas as pd
    import pyemu

    #os.chdir(os.path.join("smoother","freyberg"))

    #if not os.path.exists("freyberg.xy"):
    import flopy

    ml = flopy.modflow.Modflow.load("freyberg.nam",
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
    csv_files = [f for f in os.listdir('.') if f.endswith(".csv")]
    [os.remove(csv_file) for csv_file in csv_files]

    pst = pyemu.Pst(os.path.join("pest.pst"))
    pst.pestpp_options["ies_parameter_csv"] = "par.csv"
    pst.pestpp_options["ies_observation_csv"] = "obs.csv"
    pst.pestpp_options["ies_obs_restart_csv"] = "restart_obs.csv"
    pst.pestpp_options["parcov_filename"] = "freyberg_prior.jcb"
    #pst.observation_data.loc[pst.nnz_obs_names,"weight"] /= 10.0
    pst.write("pest.pst")

    #dia_parcov = pyemu.Cov.from_parameter_data(pst,sigma_range=6.0)
    par = pst.parameter_data
    hk_names = par.loc[par.parnme.apply(lambda x: x.startswith("hk")),'parnme']
    hk_par = par.loc[hk_names,:]
    hk_par.loc[:,"i"] = hk_par.parnme.apply(lambda x: int(x[3:5]))
    hk_par.loc[:, "j"] = hk_par.parnme.apply(lambda x: int(x[-2:]))
    hk_par.loc[:,"x"] = hk_par.apply(lambda x: ml.sr.xcentergrid[x.i,x.j],axis=1)
    hk_par.loc[:, "y"] = hk_par.apply(lambda x: ml.sr.ycentergrid[x.i, x.j], axis=1)


    gs = pyemu.utils.geostats.read_struct_file(os.path.join("structure.dat"))
    parcov = pyemu.helpers.geostatistical_prior_builder(pst,struct_dict={gs:hk_par},sigma_range=6)

    parcov.to_binary("freyberg_prior.jcb")
    parcov.to_ascii("freyberg_prior.cov")
    obscov = pyemu.Cov.from_observation_data(pst)
    num_reals = 30
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,parcov,num_reals=num_reals,use_homegrown=True)
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst,num_reals=num_reals)
    oe.to_csv("obs.csv")
    oe.to_csv("obs1.csv")
    pe.to_csv("sweep_in.csv")
    pe.to_csv("par.csv")
    pe.to_csv("par1.csv")
    pyemu.helpers.start_slaves(".","sweep","pest.pst", num_slaves=10,master_dir='.')
    df = pd.read_csv("sweep_out.csv")
    df.to_csv("restart_obs1.csv")
    df = df.loc[:,pst.observation_data.obsnme.apply(str.upper)]
    df.to_csv("restart_obs.csv")
    pst.pestpp_options["ies_parameter_csv"] = "par.csv"
    pst.pestpp_options["ies_observation_csv"] = "obs.csv"
    pst.pestpp_options["ies_obs_restart_csv"] = "restart_obs.csv"
    pst.pestpp_options["parcov_filename"] = "freyberg_prior.jcb"
    pst.write("pest.pst")

def ies():
    parcov = pyemu.Cov.from_binary("freyberg_prior.jcb")
    es = pyemu.EnsembleSmoother("pest.pst",parcov=parcov,verbose="ies.log")
    es.initialize(parensemble="par1.csv",obsensemble="obs1.csv",restart_obsensemble="restart_obs1.csv")
    es.update()
    #es.update(lambda_mults=[0.1,1.0,10.0],run_subset=10)

    #es.update(lambda_mults=[0.1,1.0,10.0],run_subset=10)

if __name__ == "__main__":
    prep()
    ies()