import os
import shutil
import numpy as np
import flopy
import pyemu

def setup():

    nlay, nrow, ncol = 3, 300, 300
    nper = 1
    perlen = 1.0
    nstp = 1
    #botm = np.arange(-10, (nlay)*10, -10)
    botm = np.linspace(-10,-100,nlay)
    m = flopy.modflow.Modflow("pest", version="mfnwt", exe_name="mfnwt", external_path='.', verbose=True)
    steady = [False for _ in range(nper)]
    # steady[0] = True
    flopy.modflow.ModflowDis(m, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper,
                             perlen=perlen, nstp=nstp, top=0.0, botm=botm, steady=steady)

    flopy.modflow.ModflowBas(m, strt=0.0, ibound=1)

    wel_step = int(nrow / 20)

    wel_i = [i for i in range(wel_step, nrow, wel_step)]
    wel_j = [j for j in range(wel_step, ncol, wel_step)]
    wel_max = 0
    wel_min = -100

    wel_data = {}
    for iper in range(nper):
        wd = []
        fluxes = np.random.uniform(wel_min, wel_max, len(wel_i))
        for i, j, f in zip(wel_i, wel_j, fluxes):
            wd.append([nlay - 1, i, j, f])
        wel_data[iper] = wd

    flopy.modflow.ModflowWel(m, stress_period_data=wel_data, ipakcb=50)

    ghb_stage1 = np.random.uniform(-1, 0, nper)
    ghb_stage2 = np.random.uniform(-1, 0, nper)
    ghb_stage_sw = np.random.uniform(-1, 0, nper)
    ghb_data = {}
    for iper, (g1, g2, gs) in enumerate(zip(ghb_stage1, ghb_stage2, ghb_stage_sw)):
        gd = []
        for i in range(0, int(nrow / 2)):
            gd.append([0, i, 0, g1, 1.0])
        for i in range(int(nrow / 2), nrow):
            gd.append([0, i, ncol - 1, g2, 10.0])
        for j in range(wel_step, ncol, wel_step):
            for i in range(nrow):
                gd.append([0, i, j, gs, 5.0])
        ghb_data[iper] = gd
    flopy.modflow.ModflowGhb(m, stress_period_data=ghb_data, ipakcb=50)

    hk = 10.0 ** (np.random.randn(nlay, nrow, ncol) * 0.01)
    vka = 10.0 ** (np.random.randn(nlay, nrow, ncol) * 0.01)
    print(hk.min(), hk.max(), hk.mean(), hk.std())
    flopy.modflow.ModflowUpw(m, hk=hk, vka=hk, ss=0.001, sy=0.1, ipakcb=50)

    flopy.modflow.ModflowOc(m, save_every=True,
                            save_types=["save head"])  # ,"save budget","print budget"])#stress_period_data=None)

    flopy.modflow.ModflowNwt(m, iprnwt=1, headtol=0.5)
    m.model_ws = "temp"
    m.write_input()
    m.run_model()
    grid_props = []
    hds_kperk = [[0,0]]
    #grid_props.append(["upw.hk",0])
    for k in range(m.nlay):
        grid_props.append(["upw.hk",k])
        grid_props.append(["upw.vka",k])
        grid_props.append(["upw.ss",k])

    ph = pyemu.helpers.PstFromFlopyModel(m,new_model_ws="template",grid_props=grid_props,hds_kperk=hds_kperk,
                                    model_exe_name="mfnwt",build_prior=False,remove_existing=True)

    ph.pst.observation_data.loc[:,"weight"] = 0.0
    ph.pst.observation_data.iloc[::5, :].loc[:,"weight"] = 1.0

    num_reals = 100
    parcov = pyemu.Cov.from_parameter_data(ph.pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(ph.pst,parcov,num_reals=num_reals)
    pe.to_csv(os.path.join("template","par.csv"))

    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(ph.pst,num_reals=num_reals)
    oe.to_csv(os.path.join("template","obs.csv"))

    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(ph.pst, num_reals=num_reals)
    oe.to_csv(os.path.join("template", "restart_obs.csv"))

def prep():
    if os.path.exists("master"):
       shutil.rmtree("master")
    shutil.copytree("template","master")
    pst = pyemu.Pst(os.path.join("master","pest.pst"))
    pst.pestpp_options["ies_parameter_csv"] = "par.csv"
    pst.pestpp_options["ies_observation_csv"] = "obs.csv"
    pst.pestpp_options["ies_obs_restart_csv"] = "restart_obs.csv"
    pst.pestpp_options["ies_use_approx"] = "true"
    pst.svd_data.eigthresh = 1.0e-5
    pst.svd_data.maxsing = 1.0e+6
    pst.control_data.noptmax = 1
    pst.write(os.path.join("master","pest.pst"))

if __name__ == "__main__":
    setup()
    prep()

