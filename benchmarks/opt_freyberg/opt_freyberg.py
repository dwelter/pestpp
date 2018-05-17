import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu

mf_nam = "freyberg.truth.nam"
mt_nam = 'freyberg.mt3d.nam'
mf_exe = "mfnwt"
mt_exe = "mt3dusgs"

new_model_ws = "template"

derinc = 2.0


def setup_models(m=None):
    org_model_ws = os.path.join("..","..","..","gw1876","models","Freyberg","Freyberg_Truth")

    if m is None:
        m = flopy.modflow.Modflow.load(mf_nam,model_ws=org_model_ws,check=False,exe_name=mf_exe)
        m.dis.nper = 1
        m.dis.perlen = 3650
        m.dis.nstp = 1
        m.dis.tsmult = 1
        m.dis.steady = True
        flopy.modflow.ModflowLmt(m,output_file_format="formatted",package_flows=["SFR"])
        m.change_model_ws("temp",reset_external=True)
        m.external_path = '.'
        m.write_input()
        m.run_model()

    mt = flopy.mt3d.Mt3dms("freyberg.mt3d",model_ws=m.model_ws,modflowmodel=m,exe_name=mt_exe,external_path='.')
    flopy.mt3d.Mt3dBtn(mt,MFStyleArr=True,prsity=0.01,sconc=0.0,icbund=m.bas6.ibound.array,perlen=3650)
    flopy.mt3d.Mt3dGcg(mt,mxiter=100)#,cclose=1.0e-7)
    #flopy.mt3d.Mt3dRct(mt,isothm=0,ireact=1,igetsc=0,rc1=0.02)
    flopy.mt3d.Mt3dAdv(mt,mixelm=-1)

    ib = m.bas6.ibound[0].array
    ssm_cells = []
    for i in range(m.nrow):
        for j in range(m.ncol):
            if ib[i,j] == 0:
                continue
            #if j >= 15: #no loading in or across the stream
            #    continue
            #ssm_cells.append([0,i,j,max(derinc,float(j+1)/10.0),15])
            ssm_cells.append([0,i,j,1.0,15])
    flopy.mt3d.Mt3dSsm(mt,crch=0.0,stress_period_data={0:ssm_cells,1:ssm_cells,2:ssm_cells})

    nstrm = np.abs(m.sfr.nstrm)
    flopy.mt3d.Mt3dSft(mt,nsfinit=nstrm,mxsfbc=0,ietsfr=0,ioutobs=1001,nobssf=nstrm,
                       obs_sf=np.arange(nstrm)+1)

    mt.write_input()
    mt.run_model()
    return mt


def setup_pest():
    m = flopy.modflow.Modflow.load(mf_nam, model_ws="temp", check=False, exe_name=mf_exe)

    props = [["upw.hk",0],["rch.rech",0],["extra.prst",0],["extra.rc11",0]]
    kperk = [[m.nper-1,0]]
    ph = pyemu.helpers.PstFromFlopyModel(mf_nam,org_model_ws="temp",new_model_ws=new_model_ws,grid_props=props,
                                         const_props=props,sfr_pars=True,remove_existing=True,
                                         model_exe_name=mf_exe,hds_kperk=kperk,
                                         extra_post_cmds=["{0} {1}".format(mt_exe,mt_nam)])

    pyemu.helpers.run("{0} {1}".format(mf_exe, mf_nam),cwd=new_model_ws)

    

    # mt = flopy.mt3d.Mt3dms.load(mt_nam, model_ws="temp", exe_name=mt_exe)
    # mt.change_model_ws(new_model_ws, reset_external=True)
    # mt.external_path = '.'
    # mt.write_input()
    # pyemu.helpers.run("{0} {1}".format(mt_exe,mt_nam),cwd=new_model_ws)
    mt = setup_models(ph.m)
    
    [shutil.copy2(os.path.join(new_model_ws, f), os.path.join(new_model_ws, 'arr_org', f)) for f in
     os.listdir(new_model_ws) if "prsity" in f.lower()]

    [shutil.copy2(os.path.join(new_model_ws, f), os.path.join(new_model_ws, 'arr_org', f)) for f in
     os.listdir(new_model_ws) if "rc11" in f.lower()]

    # mod arr pars csv for prsity name len issue - 12 chars, seriously?
    df = pd.read_csv(os.path.join(new_model_ws, "arr_pars.csv"))
    df.loc[:, "model_file"] = df.model_file.apply(lambda x: x.replace("prst", "prsity"))
    df.loc[:, "org_file"] = df.org_file.apply(lambda x: x.replace("prst", "prsity"))
    df.to_csv(os.path.join(new_model_ws,"arr_pars.csv"))

    df = write_ssm_tpl()

    os.chdir(new_model_ws)

    ph.pst.add_parameters("freyberg.mt3d.ssm.tpl","freyberg.mt3d.ssm")
    par = ph.pst.parameter_data
    for col in df.columns:
        if col in par.columns:
            par.loc[df.index,col] = df.loc[:,col]

    frun_line, ins_filenames, df = pyemu.gw_utils.setup_mtlist_budget_obs("freyberg.mt3d.list",
                                                                          start_datetime=None)
    ph.frun_post_lines.append(frun_line)
    for ins_file in ins_filenames:
        ph.pst.add_observations(ins_file,ins_file.replace(".ins",""))
    obs = ph.pst.observation_data
    for col in df.columns:
        if col in obs.columns:
            obs.loc[df.index,col] = df.loc[:,col]

    ph.tmp_files.append('mt3d_link.ftl')
    df = pyemu.gw_utils.setup_sft_obs("freyberg.mt3d",times=[np.cumsum(mt.btn.perlen.array)[-1]])
    ph.pst.add_observations("freyberg.mt3d.processed.ins","freyberg.mt3d.processed")
    ph.pst.observation_data.loc[df.index,"obgnme"] = df.index.map(lambda x: x[:4])
    ph.pst.observation_data.loc[df.index,"obsval"] = df.obsval
    ph.tmp_files.append("freyberg.mt3d")
    ph.frun_post_lines.append("pyemu.gw_utils.apply_sft_obs()")

    fline, df = pyemu.gw_utils.setup_hds_obs("MT3D001.ucn",kperk,skip=1.0e+30,prefix="ucn")
    ph.pst.add_observations("MT3D001.ucn.dat.ins","MT3D001.ucn.dat")
    ph.pst.observation_data.loc[df.index,"obgnme"] = "ucn"
    ph.pst.observation_data.loc[df.index,"obsval"] = df.obsval
    ph.tmp_files.append("MT3D001.ucn")
    ph.frun_post_lines.append(fline)

    os.chdir("..")

    ph.write_forward_run()
    ph.pst.control_data.noptmax = 0
    ph.pst.rectify_pgroups()

    ph.pst.parameter_data.sort_values(by=["pargp","parnme"],inplace=True)

    par = ph.pst.parameter_data
    #par.loc[par.pargp=="grrc110","parval1"] = 2.0

    ph.pst.parameter_groups.loc["kg","inctyp"] = "absolute"
    ph.pst.parameter_groups.loc["kg","derinc"] = derinc

    # sw conc constraint
    obs = ph.pst.observation_data
    obs.loc[:,"weight"] = 0.0
    sw_conc_obs = obs.loc[obs.obgnme=="sfrc","obsnme"]
    obs.loc[sw_conc_obs,"obgnme"] = "less_swconc"

    #only turn on one constraint in the middle of the domain
    obs.loc["sfrc30_1_03650.00","weight"] = 1.0
    obs.loc["sfrc30_1_03650.00","obsval"] *= 1.5 #% increase

    # concentration constraints at pumping wells
    wel_df = pd.DataFrame.from_records(ph.m.wel.stress_period_data[0])
    print(wel_df.dtypes)
    wel_df.loc[:,"obsnme"] = wel_df.apply(lambda x: "ucn_{0:02.0f}_{1:03.0f}_{2:03.0f}_000".format(x.k,x.i,x.j),axis=1)
    obs.loc[wel_df.obsnme,"obgnme"] = "less_wlconc"
    obs.loc[wel_df.obsnme,"weight"] = 1.0
    obs.loc[wel_df.obsnme,"obsval"] *= 1.5 #% increase


    # pumping well mass constraint
    #obs.loc["gw_we1c_003650.0","obgnme"] = "greater_well"
    #obs.loc["gw_we1c_003650.0", "weight"] = 1.0
    


    # constant head mass constraint
    # obs.loc["gw_cohe1c_003650.0","obgnme"] = "greater_ch"
    # obs.loc["gw_cohe1c_003650.0", "weight"] = 1.0
    # obs.loc["gw_cohe1c_003650.0", "obsval"] *= 1.2 #% increase

    # fix all non dec vars so we can set upper and lower bound constraints
    par = ph.pst.parameter_data
    par.loc[par.pargp!="kg","partrans"] = "fixed"

    # add some pi constraints to make sure all dec vars have at
    # least one element in the response matrix
    # set lower bounds
    parval1 = par.parval1.copy()
    par.loc[par.pargp == "kg", "parval1"] = par.loc[par.pargp == "kg", "parlbnd"]
    pyemu.helpers.zero_order_tikhonov(pst=ph.pst)
    ph.pst.prior_information.loc[:,"obgnme"] = "greater_bnd"

    # set upper bounds
    par.loc[par.pargp=="kg","parval1"] = par.loc[par.pargp=="kg","parubnd"]
    pyemu.helpers.zero_order_tikhonov(ph.pst,reset=False)
    ph.pst.prior_information.loc[:, "obgnme"] = "less_bnd"


    #set dec vars to background value (derinc)
    par.loc[par.pargp == "kg", "parval1"] = parval1

    #unfix pars
    par.loc[par.pargp != "kg", "partrans"] = "log"

    ph.pst.pestpp_options = {}
    ph.pst.pestpp_options["opt_dec_var_groups"] = "kg"
    ph.pst.pestpp_options["opt_obj_func"] = 'obj.dat'
    ph.pst.pestpp_options["opt_direction"] = "max"
    #ph.pst.pestpp_options["opt_risk"] = 0.95
    ph.pst.write(os.path.join(new_model_ws,"freyberg.pst"))
    pyemu.helpers.run("pestpp freyberg.pst",cwd=new_model_ws)

    with open(os.path.join(new_model_ws,"obj.dat"),'w') as f:
        for pname in par.loc[par.pargp=="kg","parnme"]:
            f.write("{0} {1}\n".format(pname,1.0))

def write_ssm_tpl():
    ssm = os.path.join(new_model_ws,"freyberg.mt3d.ssm")
    f_in = open(ssm,'r')
    f_tpl = open(ssm+".tpl",'w')
    f_tpl.write("ptf ~\n")
    parnme,parval1 = [],[]
    while True:
        line = f_in.readline().lower()
        if line == "":
            break
        f_tpl.write(line)
        if "stress period" in line:
            nss = int(line.strip().split()[0])
            for i in range(nss):
                line = f_in.readline().lower()
                if line == '':
                    raise Exception()
                raw = line.strip().split()
                l,r,c = [int(r) for r in raw[:3]]
                parval1.append(float(raw[3]))
                pn = "~k_{0:02d}~".format(r-1)
                #pn = "~k_{0:02d}_{1:02d}~".format(r-1,c-1)
                line = " {0:9d} {1:9d} {2:9d} {3:9s} {4:9d}\n".format(l,r,c,pn,15)
                f_tpl.write(line)
                parnme.append(pn)
            f_tpl.write("-1\n-1\n-1\n-1\n-1\n")
            break
    df = pd.DataFrame({"tpl_str":parnme})
    df.loc[:,'parnme'] = df.tpl_str.apply(lambda x: x.replace("~",""))
    df.index = df.parnme
    df.loc[:,"pargp"] = "kg"
    df.loc[:,"parval1"] = parval1
    #df.loc[:,"parubnd"] = df.parval1 * 2.0
    df.loc[:,"parubnd"] = df.parval1 * 10000000.0
    df.loc[:,"partrans"] = "none"
    df.loc[:,'parlbnd'] = df.parval1 * 0.9
    #df.loc[:,"parlbnd"] = 0.0 #here we let loading decrease...
    #df.loc[:,"parval1"] = derinc
    #df.loc[:,"j"] = df.parnme.apply(lambda x: int(x.split('_')[2]))

    return df


def run_jco():
    pst = pyemu.Pst(os.path.join(new_model_ws,"freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(new_model_ws,"freyberg.pst"))

    pyemu.helpers.start_slaves(new_model_ws,"pestpp","freyberg.pst",num_slaves=15,slave_root='.',
                               master_dir="resp_master")

def run_pestpp_opt():
    pst_file = os.path.join(new_model_ws, "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.control_data.noptmax = 1
    pst.write(pst_file)
    pyemu.helpers.start_slaves(new_model_ws,"pestpp-opt","freyberg.pst",num_slaves=15,master_dir="master_opt",
                               slave_root='.')

def spike_test():
    pst_file = os.path.join(new_model_ws, "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.parameter_data.loc["k_00","parval1"] += derinc
    #par = pst.parameter_data
    #fosm_pars = par.loc[par.pargp != "kg", "parnme"]
    #par.loc[fosm_pars,"parval1"] = par.loc[fosm_pars,"parubnd"]
    pst.control_data.noptmax = 0
    pst.write(pst_file)
    pyemu.helpers.run("pestpp {0}".format(os.path.split(pst_file)[-1]),cwd=new_model_ws)
    pst = pyemu.Pst(pst_file)
    pst.res.loc[:,"wr"] = pst.res.residual * pst.res.weight
    print(pst.res.loc[pst.res.wr>0,:])
    print(pst.res.loc[pst.nnz_obs_names,:])

    m = flopy.modflow.Modflow.load(mf_nam, model_ws=new_model_ws, check=False, exe_name=mf_exe)

    obs = pst.observation_data
    ucn_obs = obs.loc[obs.obsnme.apply(lambda x: x.startswith("ucn")),:]
    ucn_obs.loc[:,"i"] = ucn_obs.obsnme.apply(lambda x: int(x.split("_")[2]))
    ucn_obs.loc[:, "j"] = ucn_obs.obsnme.apply(lambda x: int(x.split("_")[3]))

    fig = plt.figure(figsize=(30,10))
    ax1 = plt.subplot(131)
    arr = np.zeros((m.nrow, m.ncol))
    arr[ucn_obs.i, ucn_obs.j] = pst.res.loc[ucn_obs.obsnme, "measured"]
    arr = np.ma.masked_where(m.bas6.ibound[0].array == 0, arr)
    cb = ax1.imshow(arr)
    plt.colorbar(cb)

    ax2 = plt.subplot(132)
    arr = np.zeros((m.nrow, m.ncol))
    arr[ucn_obs.i, ucn_obs.j] = pst.res.loc[ucn_obs.obsnme, "modelled"]
    arr = np.ma.masked_where(m.bas6.ibound[0].array == 0, arr)
    cb = ax2.imshow(arr)
    plt.colorbar(cb)

    ax3 = plt.subplot(133)
    arr = np.zeros((m.nrow, m.ncol))
    arr[ucn_obs.i, ucn_obs.j] = pst.res.loc[ucn_obs.obsnme, "residual"]
    arr = np.ma.masked_where(m.bas6.ibound[0].array == 0, arr)
    cb = ax3.imshow(arr)
    plt.colorbar(cb)

    plt.show()


def jco_invest():
    pst_file = os.path.join("master_opt","freyberg.pst")
    pst = pyemu.Pst(pst_file)
    par = pst.parameter_data
    fosm_pars = par.loc[par.pargp!="kg","parnme"]
    jco = pyemu.Jco.from_binary(pst_file.replace(".pst",".1.jcb")).to_dataframe()
    #print(jco.columns)
    #jco = jco.loc[:,jco.columns.map(lambda x: "kg" in x)]
    #fosm_jco = jco.loc[:,fosm_pars]
    #print(fosm_jco.loc[pst.nnz_obs_names,:])
    for oname in pst.nnz_obs_names:

        print(oname,"\n",jco.loc[oname,:])

def run_risk_sweep():
    sw_conc_inc = 1.5
    wel_mass_inc = 1.5
    chd_mass_inc = 1.5
    sw_conc_const = "sfrc20_1_03650.00"
    wel_mass_const = "gw_we1c_003650.0"
    chd_mass_const = "gw_cohe1c_003650.0"

    results_d = "results_test"
    if os.path.exists(results_d):
        shutil.rmtree(results_d)
    os.mkdir(results_d)
    restart_d = "test"
    if os.path.exists(restart_d):
        shutil.rmtree(restart_d)
    shutil.copytree(new_model_ws, restart_d)
    for f in os.listdir("_restart_files"):
        shutil.copy2(os.path.join("_restart_files", f), os.path.join(restart_d, f))

    pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = "true"

    # for inc,nme in zip([sw_conc_inc,wel_mass_inc,chd_mass_inc],
    #                    [sw_conc_const,wel_mass_const,chd_mass_const]):
    #     pst.observation_data.loc[nme,"obsval"] *= inc

    obs = pst.observation_data
    obs.loc[sw_conc_const,"obsval"] *= sw_conc_inc

    pst.control_data.noptmax = 1
    risks = np.arange(0.01, 1.0, 0.01)

    if os.path.exists("results"):
        shutil.rmtree("results")
    os.mkdir("results")

    phis, infeas = [], []
    for risk in risks:
        pst.pestpp_options["opt_risk"] = risk
        pst.write(os.path.join(restart_d, "freyberg.pst"))
        pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=restart_d)
        inf, phi = scrape_recfile(os.path.join(restart_d, "freyberg.rec"))
        phis.append(phi)
        infeas.append(inf)

    df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
    df.to_csv(os.path.join(results_d, "base.csv"))
    df.loc[df.infeas,"phi"] = np.NaN
    ax = plt.subplot(111)
    ax.plot(df.risk,df.phi)
    ax.set_xlabel("risk")
    ax.set_ylabel("$\phi$")
    ax.grid()
    plt.savefig("risk.pdf")

def run_risk_sweep_pargp():
    sw_conc_inc = 1.1
    wel_mass_inc = 1.1
    chd_mass_inc = 1.1
    results_d = "results_test"
    restart_d = "test"
    if os.path.exists(restart_d):
        shutil.rmtree(restart_d)
    shutil.copytree(new_model_ws, restart_d)
    for f in os.listdir("_restart_files"):
        shutil.copy2(os.path.join("_restart_files",f),os.path.join(restart_d,f))

    pst = pyemu.Pst(os.path.join(new_model_ws,"freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = "true"

    pst.control_data.noptmax = 1
    risks = np.arange(0.01,1.0,0.01)

    if os.path.exists("results"):
        shutil.rmtree("results")
    os.mkdir("results")

    phis, infeas = [], []
    for risk in risks:
        pst.pestpp_options["opt_risk"] = risk
        pst.write(os.path.join(restart_d, "freyberg.pst"))
        pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=restart_d)
        inf, phi = scrape_recfile(os.path.join(restart_d, "freyberg.rec"))
        phis.append(phi)
        infeas.append(inf)

    df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
    df.to_csv(os.path.join("results", "base.csv"))

    jco = pyemu.Jco.from_binary(os.path.join(restart_d,"restart.jcb"))
    for pargp in pst.par_groups:
        if pargp == "kg":
            continue
        pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
        pst.pestpp_options["base_jacobian"] = "restart_pargp.jcb"
        pst.pestpp_options["hotstart_resfile"] = "restart.rei"
        pst.pestpp_options["opt_skip_final"] = "true"
        pst.control_data.noptmax = 1
        par = pst.parameter_data
        par.loc[par.pargp==pargp,"partrans"] = "fixed"
        jco_pargp = jco.get(col_names=pst.adj_par_names)
        jco_pargp.to_binary(os.path.join(restart_d,"restart_pargp.jcb"))

        phis,infeas = [],[]
        for risk in risks:
            pst.pestpp_options["opt_risk"] = risk
            pst.write(os.path.join(restart_d,"freyberg.pst"))
            pyemu.helpers.run("pestpp-opt freyberg.pst",cwd=restart_d)
            inf,phi = scrape_recfile(os.path.join(restart_d,"freyberg.rec"))
            phis.append(phi)
            infeas.append(inf)



        df = pd.DataFrame({"risk":risks,"phi":phis,"infeas":infeas})
        df.to_csv(os.path.join("results","{0}.csv".format(pargp)))
        #break

    # df.loc[df.infeas,"phi"] = np.NaN
    # ax = plt.subplot(111)
    # ax.plot(df.risk,df.phi)
    # ax.set_xlabel("risk")
    # ax.set_ylabel("$\phi$")
    # ax.grid()
    # plt.savefig("risk.pdf")

def run_risk_sweep_obgnme():
    restart_d = "test"
    if os.path.exists(restart_d):
        shutil.rmtree(restart_d)
    shutil.copytree(new_model_ws, restart_d)
    for f in os.listdir("_restart_files"):
        shutil.copy2(os.path.join("_restart_files", f), os.path.join(restart_d, f))

    pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = "true"

    pst.control_data.noptmax = 1
    risks = np.arange(0.01, 1.0, 0.01)

    if os.path.exists("results_obs"):
        shutil.rmtree("results_obs")
    os.mkdir("results_obs")

    phis, infeas = [], []
    for risk in risks:
        pst.pestpp_options["opt_risk"] = risk
        pst.write(os.path.join(restart_d, "freyberg.pst"))
        pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=restart_d)
        inf, phi = scrape_recfile(os.path.join(restart_d, "freyberg.rec"))
        phis.append(phi)
        infeas.append(inf)

    df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
    df.to_csv(os.path.join("results_obs", "base.csv"))

    jco = pyemu.Jco.from_binary(os.path.join(restart_d, "restart.jcb"))
    for gp in pst.obs_groups:
        if gp.startswith("less") or gp.startswith("great"):
            continue
        pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
        pst.pestpp_options["base_jacobian"] = "restart.jcb"
        pst.pestpp_options["hotstart_resfile"] = "restart.rei"
        pst.pestpp_options["opt_skip_final"] = "true"
        pst.control_data.noptmax = 1
        obs = pst.observation_data
        obs.loc[obs.obgnme==gp,"weight"] = 1.0

        phis, infeas = [], []
        for risk in risks:
            pst.pestpp_options["opt_risk"] = risk
            pst.write(os.path.join(restart_d, "freyberg.pst"))
            pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=restart_d)
            inf, phi = scrape_recfile(os.path.join(restart_d, "freyberg.rec"))
            phis.append(phi)
            infeas.append(inf)

        df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
        df.to_csv(os.path.join("results_obs", "{0}.csv".format(pargp)))

def scrape_recfile(recfile):
    infeas = False
    with open(recfile,'r') as f:
        for line in f:
            if "best objective function value" in line:
                phi = float(line.strip().split(':')[-1])
                #break

            if "warning: primal solution infeasible" in line:
                infeas = True
    return infeas, phi

def plot_risk_sweep():
    dfs = {f.split('.')[0]:pd.read_csv(os.path.join("results",f)) for f in os.listdir("results")}
    ax = plt.subplot(111)
    for pargp,df in dfs.items():
        df.loc[df.infeas,"phi"] = np.NaN
        ax.plot(df.risk,df.phi,label=pargp)
    ax.grid()
    ax.legend()
    plt.savefig("risk_sweep.pdf")

def plot_loading(parfile=None):
    if parfile is not None:
        df = pyemu.pst_utils.read_parfile(parfile)
    else:
        df = pyemu.pst_utils.read_parfile(os.path.join("master_opt","freyberg.par"))

    df = df.loc[df.parnme.apply(lambda x: x.startswith("k_")),:]
    print(df)


if __name__ == "__main__":
    #write_ssm_tpl()
    #run_test()
    setup_models()
    setup_pest()
    #spike_test()

    run_pestpp_opt()
    #jco_invest()
    #run_risk_sweep()
    #plot_loading()
    #plot_risk_sweep()
    #run_risk_sweep_obgnme()
    #run_risk_sweep()

