import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import flopy
import pyemu

mf_nam = "freyberg.truth.nam"
mt_nam = 'freyberg.mt3d.nam'
mf_exe = "mfnwt"
mt_exe = "mt3dusgs"

new_model_ws = "template"

derinc = 1.0

phi_factor = 365. * 40.

org_model_ws = os.path.join("..","..","..","gw1876","models","Freyberg","Freyberg_Truth")


def setup_truth():
    setup_models()
    if os.path.exists("truth"):
        shutil.rmtree("truth")
    shutil.copytree("temp","truth")
    shutil.copy2(os.path.join(org_model_ws,"hk.zones"),os.path.join("truth","hk.zones"))
    zarr = np.loadtxt(os.path.join("truth","hk.zones"),dtype=np.int)
    uvals = np.unique(zarr)
    print(uvals)
    rc1_range = (0.0001,0.01)
    prsity_range = (0.01,0.2)
    hk_range = (1.0,100.0)
    ss_range = (0.00001,0.001)
    sy_range = (0.01,0.1)
    #rch_range


    

def setup_models(m=None):
    
    if m is None:
        m = flopy.modflow.Modflow.load(mf_nam,model_ws=org_model_ws,check=False,exe_name=mf_exe)
        m.dis.nper = 3
        m.dis.perlen = [3650,3650,3650]
        m.dis.nstp = 1
        m.dis.tsmult = 1
        m.dis.steady = [True,False,False]
        m.sfr.reach_data['j'] += 1
        flopy.modflow.ModflowLmt(m,output_file_format="formatted",package_flows=["SFR"])
        m.change_model_ws("temp",reset_external=True)
        m.external_path = '.'
        m.write_input()
        m.run_model()
        for f in ["obs_loc.csv","hk.zones"]:
            shutil.copy2(os.path.join(org_model_ws,f),os.path.join("temp",f))

    mt = flopy.mt3d.Mt3dms("freyberg.mt3d",model_ws=m.model_ws,modflowmodel=m,exe_name=mt_exe,external_path='.')
    flopy.mt3d.Mt3dBtn(mt,MFStyleArr=True,prsity=0.01,sconc=1.0,icbund=m.bas6.ibound.array,perlen=3650)
    flopy.mt3d.Mt3dGcg(mt,mxiter=100)#,cclose=1.0e-7)
    flopy.mt3d.Mt3dRct(mt,isothm=0,ireact=1,igetsc=0,rc1=0.001)
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
            ssm_cells.append([0,i,j,derinc,15])
    flopy.mt3d.Mt3dSsm(mt,crch=0.0,stress_period_data={0:ssm_cells,1:ssm_cells,2:ssm_cells})

    nstrm = np.abs(m.sfr.nstrm)
    flopy.mt3d.Mt3dSft(mt,nsfinit=nstrm,mxsfbc=0,ietsfr=0,ioutobs=1001,nobssf=nstrm,
                       obs_sf=np.arange(nstrm)+1)

    mt.write_input()
    mt.run_model()
    return mt


def setup_pest():
    m = flopy.modflow.Modflow.load(mf_nam, model_ws="temp", check=False, exe_name=mf_exe)
    props = [["upw.hk",0],["rch.rech",0],["rch.rech",1],["rch.rech",2],["extra.prst",0],["extra.rc11",0],["extra.sc1",0]]
    kperk = [[0,0],[1,0],[2,0]]
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

    [shutil.copy2(os.path.join(new_model_ws, f), os.path.join(new_model_ws, 'arr_org', f)) for f in
     os.listdir(new_model_ws) if "sconc1" in f.lower()]

    # mod arr pars csv for prsity name len issue - 12 chars, seriously?
    df = pd.read_csv(os.path.join(new_model_ws, "arr_pars.csv"))
    df.loc[:, "model_file"] = df.model_file.apply(lambda x: x.replace("prst", "prsity"))
    df.loc[:, "org_file"] = df.org_file.apply(lambda x: x.replace("prst", "prsity"))

    df.loc[:, "model_file"] = df.model_file.apply(lambda x: x.replace("sc1", "sconc1"))
    df.loc[:, "org_file"] = df.org_file.apply(lambda x: x.replace("sc1", "sconc1"))
    
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
    df = pyemu.gw_utils.setup_sft_obs("freyberg.mt3d",times=np.cumsum(mt.btn.perlen.array))
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
    ph.pst.pestpp_options["parcov"] = "freyberg.truth.pst.prior.cov"
    ph.pst.control_data.noptmax = 0
    ph.pst.rectify_pgroups()

    ph.pst.parameter_data.sort_values(by=["pargp","parnme"],inplace=True)

    par = ph.pst.parameter_data
    #par.loc[par.pargp=="grrc110","parval1"] = 2.0

    ph.pst.parameter_groups.loc["kg","inctyp"] = "absolute"
    ph.pst.parameter_groups.loc["kg","derinc"] = derinc

    # hds obs - end of first sp
    obs_df = pd.read_csv(os.path.join("temp","obs_loc.csv"))
    obs_df.loc[:,"obsnme"] = obs_df.apply(lambda x: "hds_00_{0:03d}_{1:03d}_000".format(x.row-1,x.col-1),axis=1)
    print(obs_df)
    ph.pst.observation_data.loc[obs_df.obsnme,"weight"] = 1.0
    ph.pst.observation_data.loc[obs_df.obsnme,"ogbnme"] = "cal_hds"

    # sw conc constraint
    obs = ph.pst.observation_data
    obs.loc[:,"weight"] = 0.0
    sw_conc_obs = obs.loc[obs.obgnme=="sfrc","obsnme"]
    obs.loc[sw_conc_obs,"obgnme"] = "less_swconc"

    #only turn on one constraint in the middle of the domain - end of 2nd sp
    obs.loc["sfrc30_1_07300.00","weight"] = 1.0
    obs.loc["sfrc30_1_07300.00","obsval"] *= 1.2 #% increase

    # concentration constraints at pumping wells - end of 2nd sp
    wel_df = pd.DataFrame.from_records(ph.m.wel.stress_period_data[0])
    print(wel_df.dtypes)
    wel_df.loc[:,"obsnme"] = wel_df.apply(lambda x: "ucn_{0:02.0f}_{1:03.0f}_{2:03.0f}_001".format(x.k,x.i,x.j),axis=1)
    obs.loc[wel_df.obsnme,"obgnme"] = "less_wlconc"
    obs.loc[wel_df.obsnme,"weight"] = 1.0
    obs.loc[wel_df.obsnme,"obsval"] *= 1.2 #% increase

    # pumping well mass constraint
    #obs.loc["gw_we1c_003650.0","obgnme"] = "greater_well"
    #obs.loc["gw_we1c_003650.0", "weight"] = 1.0

    # constant head mass constraint
    # obs.loc["gw_cohe1c_003650.0","obgnme"] = "greater_ch"
    # obs.loc["gw_cohe1c_003650.0", "weight"] = 1.0
    # obs.loc["gw_cohe1c_003650.0", "obsval"] *= 1.2 #% increase

    # pestpp-opt does this internally now...
    # fix all non dec vars so we can set upper and lower bound constraints
    # par = ph.pst.parameter_data
    # par.loc[par.pargp!="kg","partrans"] = "fixed"

    # # add some pi constraints to make sure all dec vars have at
    # # least one element in the response matrix
    # # set lower bounds
    # parval1 = par.parval1.copy()
    # par.loc[par.pargp == "kg", "parval1"] = par.loc[par.pargp == "kg", "parlbnd"]
    # pyemu.helpers.zero_order_tikhonov(pst=ph.pst)

    # ph.pst.prior_information.loc[:,"pilbl"] += "_l"
    # l_const = ph.pst.prior_information.pilbl.copy()
    # #ph.pst.prior_information.loc[:,"obgnme"] = "greater_bnd"

    # # set upper bounds
    # par.loc[par.pargp=="kg","parval1"] = par.loc[par.pargp=="kg","parubnd"]
    # pyemu.helpers.zero_order_tikhonov(ph.pst,reset=False)
    # ph.pst.prior_information.loc[:, "obgnme"] = "less_bnd"
    # ph.pst.prior_information.loc[l_const,"obgnme"] = "greater_bnd"


    #set dec vars to background value (derinc)
    #par.loc[par.pargp == "kg", "parval1"] = parval1

    #unfix pars
    #par.loc[par.pargp != "kg", "partrans"] = "log"

    ph.pst.pestpp_options = {}
    ph.pst.pestpp_options["opt_dec_var_groups"] = "k1"
    ph.pst.pestpp_options["opt_obj_func"] = 'obj.dat'
    ph.pst.pestpp_options["opt_direction"] = "max"
    #ph.pst.pestpp_options["opt_risk"] = 0.95
    ph.pst.write(os.path.join(new_model_ws,"freyberg.pst"))
    pyemu.helpers.run("pestpp freyberg.pst",cwd=new_model_ws)

    with open(os.path.join(new_model_ws,"obj.dat"),'w') as f:
        for pname in par.loc[par.pargp=="k1","parnme"]:
            f.write("{0} {1}\n".format(pname,1.0))

def write_ssm_tpl():
    ssm = os.path.join(new_model_ws,"freyberg.mt3d.ssm")
    f_in = open(ssm,'r')
    f_tpl = open(ssm+".tpl",'w')
    f_tpl.write("ptf ~\n")
    parnme,parval1 = [],[]
    kper = -1
    while True:
        line = f_in.readline().lower()
        if line == "":
            break
        f_tpl.write(line)
        if "stress period" in line:
            kper += 1
            nss = int(line.strip().split()[0])
            for i in range(nss):
                line = f_in.readline().lower()
                if line == '':
                    raise Exception()
                raw = line.strip().split()
                l,r,c = [int(r) for r in raw[:3]]
                parval1.append(float(raw[3]))
                #pn = "~k_all~"
                #pn = "~k_{0:02d}~".format(r-1)
                pn = "~k{0:1d}_{1:02d}_{2:02d}~".format(kper,r-1,c-1)
                line = " {0:9d} {1:9d} {2:9d}{3:10s} {4:9d}\n".format(l,r,c,pn,15)
                f_tpl.write(line)
                parnme.append(pn)
            #f_tpl.write("-1\n-1\n-1\n-1\n-1\n")
    f_tpl.close()
    f_in.close()
    df = pd.DataFrame({"tpl_str":parnme})
    df.loc[:,'parnme'] = df.tpl_str.apply(lambda x: x.replace("~",""))
    df.index = df.parnme
    df.loc[:,"pargp"] = df.parnme.apply(lambda x: x.split('_')[0])
    df.loc[:,"parval1"] = parval1
    df.loc[:,"parubnd"] = 100.0
    #df.loc[:,"parubnd"] = df.parval1 * 10000000.0
    df.loc[:,"partrans"] = "none"
    #df.loc[:,'parlbnd'] = df.parval1 * 0.9
    df.loc[:,"parlbnd"] = 0.0 #here we let loading decrease...
    #df.loc[:,"parval1"] = derinc
    #df.loc[:,"j"] = df.parnme.apply(lambda x: int(x.split('_')[2]))

    return df

def plot_risk_sweep_pargp():
    r_d = "results_pargp1"
    files = os.listdir(r_d)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    colors=['g','m','c','b','r','y',"maroon","indigo","brown",'darkred',"olivedrab","coral"]
    for i,f in enumerate(files):
        if "base" in f:
            c = 'k'
            lw = 3
            ls = '--'
        else:
            c = colors[i]
            lw = 1.0
            ls = '-'
        df = pd.read_csv(os.path.join(r_d,f),index_col=0)
        df.loc[df.infeas==True,"phi"] = np.NaN
        df.loc[:,"phi"] *= phi_factor
        ax.plot(df.risk,df.phi.values,color=c,lw=lw,ls=ls,alpha=0.5,label=f.split('.')[0])
    ax.legend()
    ax.grid()
    ax.set_ylabel("$\\phi (\\$)$")
    ax.set_xlabel("risk")
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1e'))
    plt.savefig("risk_pargp.pdf")
    plt.show()


def plot_risk_sweep_obgnme():
    r_d = "results_obgnme"
    files = os.listdir(r_d)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    colors=['g','m','c','b','r','y',"maroon","indigo","brown",'darkred',"olivedrab","coral"]
    for i,f in enumerate(files):
        if "base" in f:
            c = 'k'
            lw = 3
            ls = '--'
        else:
            c = "0.5"
            lw = 1.0
            ls = '-'
        df = pd.read_csv(os.path.join(r_d,f),index_col=0)
        df.loc[df.infeas==True,"phi"] = np.NaN
        df.loc[:,"phi"] *= phi_factor
        ax.plot(df.risk,df.phi.values,color=c,lw=lw,ls=ls,alpha=0.5,label=f.split('.')[0])
    ax.legend()
    ax.grid()
    ax.set_ylabel("$\\phi (\\$)$")
    ax.set_xlabel("risk")
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1e'))
    plt.savefig("risk_obgnme.pdf")
    plt.show()


def run_jco():
    pst = pyemu.Pst(os.path.join(new_model_ws,"freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(new_model_ws,"freyberg.pst"))

    pyemu.helpers.start_slaves(new_model_ws,"pestpp","freyberg.pst",num_slaves=15,slave_root='.',
                               master_dir="resp_master")


def start_slaves():
    pyemu.helpers.start_slaves(new_model_ws,"pestpp-opt","freyberg.pst",num_slaves=15,master_dir=None,
                               slave_root='.',port=4005)

def run_pestpp_opt():
    pst_file = os.path.join(new_model_ws, "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.control_data.noptmax = 1
    pst.write(pst_file)
    pyemu.helpers.start_slaves(new_model_ws,"pestpp-opt","freyberg.pst",num_slaves=15,master_dir="master",
                               slave_root='.',port=4005)

def run_pestpp_gsa():
    pst_file = os.path.join(new_model_ws, "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.control_data.noptmax = 1
    par = pst.parameter_data
    par.loc[par.pargp!="kg","partrans"] = "fixed"
    print(pst.npar_adj)
    pst.write(os.path.join(new_model_ws,"freyberg_gsa.pst"))
    pyemu.helpers.start_slaves(new_model_ws,"pestpp-gsa","freyberg_gsa.pst",num_slaves=15,master_dir="master",
                               slave_root='.',port=4005)


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
    pst_file = os.path.join("_restart_files","restart.pst")
    pst = pyemu.Pst(pst_file)
    par = pst.parameter_data
    fosm_pars = par.loc[par.pargp!="kg","parnme"]
    jco = pyemu.Jco.from_binary(pst_file.replace(".pst",".jcb")).to_dataframe()
    rc_pars = par.loc[par.pargp.apply(lambda x: "rc1" in x),"parnme"]
    print(rc_pars)
    print(jco.loc[:,rc_pars])
    #print(jco.columns)
    #jco = jco.loc[:,jco.columns.map(lambda x: "kg" in x)]
    #fosm_jco = jco.loc[:,fosm_pars]
    #print(fosm_jco.loc[pst.nnz_obs_names,:])

def prep_for_risk_sweep():
    r_d = "_restart_files"
    if os.path.exists(r_d):
        shutil.rmtree(r_d)
    os.mkdir(r_d)
    d = "test"
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(new_model_ws,d)
    
    pst_file = os.path.join(new_model_ws, "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.control_data.noptmax = 0
    pst.write(os.path.join(d,"restart.pst"))
    pyemu.os_utils.run("pestpp restart.pst",cwd=d)
    shutil.copy2(os.path.join(d,"restart.rei"),os.path.join(r_d,"restart.rei"))
    pst.control_data.noptmax = 1
    pst.pestpp_options["opt_risk"] = 0.1
    pst.write(os.path.join(new_model_ws,"restart.pst"))
    pyemu.helpers.start_slaves(new_model_ws,"pestpp-opt","restart.pst",num_slaves=20,master_dir=d,
                               slave_root='.',port=4005)
    shutil.copy2(os.path.join(d,"restart.1.jcb"),os.path.join(r_d,"restart.jcb"))


def run_risk_sweep():

    results_d = "sweep_results"
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

    pst.control_data.noptmax = 1
    risks = np.arange(0.01, 1.01, 0.01)

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
    df.loc[:,"phi"] *= phi_factor
    df.to_csv(os.path.join(results_d, "base.csv"))
    df.loc[df.infeas,"phi"] = np.NaN
    ax = plt.subplot(111)
    ax.plot(df.risk,df.phi)
    ax.set_xlabel("risk")
    ax.set_ylabel("$\phi (\\$) $")
    ax.grid()
    plt.savefig("risk.pdf")


def run_risk_sweep_pargp():

    results_d = "results_pargp1"
    if os.path.exists(results_d):
        shutil.rmtree(results_d)
    os.mkdir(results_d)
    
    pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = "true"
    par = pst.parameter_data
    cn_pars = par.loc[par.pargp=="cn","parnme"]
    par.loc[cn_pars,"pargp"] = cn_pars.apply(lambda x:x.split('_')[0]+"_mean")
    par_base = pst.parameter_data.copy()
    pst.control_data.noptmax = 1
    risks = np.arange(0.1, 1.1, 0.1)
    
    def run(name, pst):
        d = "test"
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(new_model_ws, d)
        for f in os.listdir("_restart_files"):
            shutil.copy2(os.path.join("_restart_files", f), os.path.join(d, f))
        phis, infeas = [], []
        for risk in risks:
            pst.pestpp_options["opt_risk"] = risk
            pst.write(os.path.join(d, "freyberg.pst"))
            pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=d)
            inf, phi = scrape_recfile(os.path.join(d, "freyberg.rec"))
            phis.append(phi)
            infeas.append(inf)

        df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
        df.to_csv(os.path.join(results_d, "{0}.csv".format(name)))

    run("base",pst)

    # for pargp in pst.par_groups:
    #     if pargp == "kg":
    #         continue
    #     par = par_base.copy()
    #     par.loc[par.pargp==pargp,"partrans"] = "fixed"
    #     pst.parameter_data = par
    #     run(pargp,pst)
    for prefix in ["hk","rc1","grrech0","prst","sc1","k0","grrech1","hcond1"]:
        par = par_base.copy()
        par.loc[par.pargp.apply(lambda x:prefix in x),"partrans"] = "fixed"

        pst.parameter_data = par
        #print(pst.parameter_data.loc[pst.parameter_data.partrans == "fixed", :])
        run(prefix,pst)




def run_risk_sweep_obgnme():

    results_d = "results_obgnme"
    if os.path.exists(results_d):
        shutil.rmtree(results_d)
    os.mkdir(results_d)
    
    pst = pyemu.Pst(os.path.join(new_model_ws, "freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = "true"
    pst.control_data.noptmax = 1
    risks = np.arange(0.1, 1.1, 0.1)
    
    def run(name, pst):
        d = "test"
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(new_model_ws, d)
        for f in os.listdir("_restart_files"):
            shutil.copy2(os.path.join("_restart_files", f), os.path.join(d, f))
        phis, infeas = [], []
        for risk in risks:
            pst.pestpp_options["opt_risk"] = risk
            pst.write(os.path.join(d, "freyberg.pst"))
            pyemu.helpers.run("pestpp-opt freyberg.pst", cwd=d)
            inf, phi = scrape_recfile(os.path.join(d, "freyberg.rec"))
            phis.append(phi)
            infeas.append(inf)

        df = pd.DataFrame({"risk": risks, "phi": phis, "infeas": infeas})
        df.to_csv(os.path.join(results_d, "{0}.csv".format(name)))

    run("base",pst)

    # for pargp in pst.par_groups:
    #     if pargp == "kg":
    #         continue
    #     par = par_base.copy()
    #     par.loc[par.pargp==pargp,"partrans"] = "fixed"
    #     pst.parameter_data = par
    #     run(pargp,pst)
    obs_base = pst.observation_data.copy()
    for obgnme in pst.obs_groups:
        if "less_" in obgnme:
            continue
        obs = obs_base.copy()
        obs.loc[obs.obgnme==obgnme,"weight"] = 1.0

        pst.observation_data = obs
        #print(pst.parameter_data.loc[pst.parameter_data.partrans == "fixed", :])
        run(obgnme,pst)



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

def plot_loading(d="master"):
    pst = pyemu.Pst(os.path.join("template","freyberg.pst"))
    m = flopy.modflow.Modflow.load("freyberg.truth.nam",model_ws="template",load_only=[],check=False)
    print(pst.nnz_obs_names)
    nz_obs = pst.observation_data.loc[pst.nnz_obs_names,:]
    nz_obs = nz_obs.loc[nz_obs.obsnme.apply(lambda x: x.startswith("ucn")),:]
    nz_obs.loc[:,"i"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    nz_obs.loc[:,"j"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    nz_obs.loc[:,"x"] = nz_obs.apply(lambda x: m.sr.xcentergrid[x.i,x.j],axis=1)

    nz_obs.loc[:,"y"] = nz_obs.apply(lambda x: m.sr.ycentergrid[x.i,x.j],axis=1)
    seg_i = int([o for o in pst.nnz_obs_names if o.startswith("sfr")][0].split('_')[0][-2:])
    print(nz_obs)
    #print(seg_j)
    seg_j = 16
    seg_x,seg_y = m.sr.xcentergrid[seg_i,seg_j],m.sr.ycentergrid[seg_i,seg_j]
    dfs = [] 
    mx,mn = -1.0e+10,1.0e+10
    for i in range(pst.control_data.noptmax):
        par_file = os.path.join(d,"freyberg.{0}.par".format(i+1))
        if not os.path.exists(par_file):
            break
        df = pyemu.pst_utils.read_parfile(par_file)

        df = df.loc[df.parnme.apply(lambda x: x.startswith("k1_")),:]
        #print(df)
        df.loc[:,"i"] = df.parnme.apply(lambda x: int(x.split('_')[1]))
        df.loc[:,"j"] = df.parnme.apply(lambda x: int(x.split('_')[2]))
        df.loc[df.parval1<0,"parval1"] = 0.0
        mx = max(mx,df.parval1.max())
        mn = min(mn,df.parval1.min())
        dfs.append(df)
    for df in dfs:
        arr = np.zeros((40,20)) - 1
        arr[df.i,df.j] = df.parval1
        arr = np.ma.masked_where(arr<0,arr)
        ax = plt.subplot(111)

        cb = ax.imshow(arr,extent=m.sr.get_extent(),vmin=mn,vmax=mx)
        c = plt.colorbar(cb)
        c.set_label("N loading")
        ax.scatter(nz_obs.x,nz_obs.y,marker='x',color='r')
        ax.scatter([seg_x],[seg_y],marker='^',color='r')
        plt.show()


def constraint_testing():
    d = "test"
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree("template",d)
    for f in os.listdir("_restart_files"):
        shutil.copy2(os.path.join("_restart_files",f),os.path.join(d,f))
    m = flopy.modflow.Modflow.load("freyberg.truth.nam",model_ws="template",load_only=["wel"],check=False)

    pst = pyemu.Pst(os.path.join("template","freyberg.pst"))
    par = pst.parameter_data
    par.loc[par.pargp=="kg","parubnd"] = 10.0
    par.loc[par.pargp=="kg","parlbnd"] = 0.0
    par.loc[par.pargp=="kg","parval1"] = 5.0
    par.loc[par.pargp=="sc1","parval1"] = 5.0
    par.loc[par.pargp=="sc1","parlbnd"] = 1.0
    par.loc[par.pargp=="sc1","parubnd"] = 10.0


    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = True
    # sw conc constraint
    obs = pst.observation_data
    obs.loc[:,"weight"] = 0.0
    sw_conc_obs = obs.loc[obs.obgnme=="sfrc","obsnme"]
    obs.loc[sw_conc_obs,"obgnme"] = "less_swconc"

    #only turn on one constraint in the middle of the domain
    obs.loc["sfrc30_1_03650.00","weight"] = 1.0
    obs.loc["sfrc30_1_03650.00","obsval"] *= 1.0 #% increase

    #total SW N mass leaving the model
    obs.loc["sw_st1c_003650.0","weight"] = 1.0
    obs.loc["sw_st1c_003650.0","obsval"] *= 1.0
    obs.loc["sw_st1c_003650.0","obgnme"] = "less_mass"

    # concentration constraints at pumping wells
    wel_df = pd.DataFrame.from_records(m.wel.stress_period_data[0])
    print(wel_df.dtypes)
    wel_df.loc[:,"obsnme"] = wel_df.apply(lambda x: "ucn_{0:02.0f}_{1:03.0f}_{2:03.0f}_000".format(x.k,x.i,x.j),axis=1)
    obs.loc[wel_df.obsnme,"obgnme"] = "less_wlconc"
    obs.loc[wel_df.obsnme,"weight"] = 1.0
    obs.loc[wel_df.obsnme,"obsval"] *= 1.0 #% increase


    pst.write(os.path.join(d,"freyberg.pst"))
    pyemu.os_utils.run("pestpp-opt freyberg.pst",cwd=d)
    plot_loading(d=d)


def bounds_testing():
    m = flopy.modflow.Modflow.load("freyberg.truth.nam",model_ws="template",load_only=[],check=False)
    
    files = os.listdir("_restart_files")
    for f in files:
        shutil.copy2(os.path.join("_restart_files",f),os.path.join("master",f))
    pst = pyemu.Pst(os.path.join("master","freyberg.pst"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["opt_skip_final"] = True
    nz_obs = pst.observation_data.loc[pst.nnz_obs_names,:]
    nz_obs = nz_obs.loc[nz_obs.obsnme.apply(lambda x: x.startswith("ucn")),:]
    nz_obs.loc[:,"i"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    nz_obs.loc[:,"j"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    nz_obs.loc[:,"x"] = nz_obs.apply(lambda x: m.sr.xcentergrid[x.i,x.j],axis=1)

    nz_obs.loc[:,"y"] = nz_obs.apply(lambda x: m.sr.ycentergrid[x.i,x.j],axis=1)
    seg_i = int([o for o in pst.nnz_obs_names if o.startswith("sfr")][0].split('_')[0][-2:])
    print(nz_obs)
    #print(seg_j)
    seg_j = 16
    seg_x,seg_y = m.sr.xcentergrid[seg_i,seg_j],m.sr.ycentergrid[seg_i,seg_j]
    par_file = os.path.join("master","freyberg.1.par")
    
    print(pst.nnz_obs_names)
    mn = derinc
    mx = derinc * 1.5
    infs,phis = [],[]
    for i,ub in enumerate(np.linspace(mn,mx,100)):
        
        print(pst.par_groups)
        par = pst.parameter_data
        par.loc[par.pargp=="kg","parubnd"] = ub
        par.loc[par.pargp=="kg","parlbnd"] = 1.0

        pst.write(os.path.join("master","freyberg.pst"))
        pyemu.os_utils.run("pestpp-opt freyberg.pst",cwd="master")
        inf, phi = scrape_recfile(os.path.join("master", "freyberg.rec"))
        infs.append(inf)
        phis.append(phi)  
        df = pyemu.pst_utils.read_parfile(par_file)

        df = df.loc[df.parnme.apply(lambda x: x.startswith("k_")),:]
        #print(df)
        df.loc[:,"i"] = df.parnme.apply(lambda x: int(x.split('_')[1]))
        df.loc[:,"j"] = df.parnme.apply(lambda x: int(x.split('_')[2]))
        df.loc[df.parval1<0,"parval1"] = 0.0
        #mx = max(mx,df.parval1.max())
        #mn = min(mn,df.parval1.min())
        arr = np.zeros((40,20)) - 1
        arr[df.i,df.j] = df.parval1
        arr = np.ma.masked_where(arr<0,arr)
        ax = plt.subplot(111)

        cb = ax.imshow(arr,extent=m.sr.get_extent())#,vmin=mn,vmax=mx)
        c = plt.colorbar(cb)
        c.set_label("N loading")
        ax.scatter(nz_obs.x,nz_obs.y,marker='x',color='r')
        ax.scatter([seg_x],[seg_y],marker='^',color='r')
        ax.set_title("phi:{0:15.3G}, ub:{1:15.3G}".format(phi,ub))
        plt.savefig(os.path.join("master","load_{0:03d}.png".format(i)))
        plt.close("all")
    df = pd.DataFrame({"phi":phis,"inf":infs})
    #df.loc[df.inf==False,:] = np.NaN
    df.phi.plot()
    plt.show()

if __name__ == "__main__":
    #setup_truth()
    #write_ssm_tpl()
    #run_test()
    setup_models()
    setup_pest()
    #spike_test()
    #start_slaves()
    #run_pestpp_opt()
    
    prep_for_risk_sweep()
    #jco_invest()
    #run_risk_sweep()
    #run_risk_sweep_pargp()
    #plot_risk_sweep_pargp()
    #bounds_test()
    #plot_loading()
    #plot_risk_sweep()
    #run_risk_sweep_obgnme()
    #plot_risk_sweep_obgnme()
    #run_risk_sweep()
    #constraint_testing()

