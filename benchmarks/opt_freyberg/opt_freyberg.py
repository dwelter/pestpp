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

derinc = 0.1

def setup_models():
    org_model_ws = os.path.join("..","..","..","gw1876","models","Freyberg","Freyberg_Truth")


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

    mt = flopy.mt3d.Mt3dms("freyberg.mt3d",model_ws="temp",modflowmodel=m,exe_name=mt_exe,external_path='.')
    flopy.mt3d.Mt3dBtn(mt,MFStyleArr=True,prsity=0.01,sconc=0.0,icbund=m.bas6.ibound.array,perlen=3650)
    flopy.mt3d.Mt3dGcg(mt)
    flopy.mt3d.Mt3dRct(mt,isothm=0,ireact=0,igetsc=0,rc1=0.02)
    flopy.mt3d.Mt3dAdv(mt,mixelm=-1)

    ib = m.bas6.ibound[0].array
    ssm_cells = []
    for i in range(m.nrow):
        for j in range(m.ncol):
            if ib[i,j] == 0:
                continue
            ssm_cells.append([0,i,j,0.0,15])
    flopy.mt3d.Mt3dSsm(mt,crch=0.0,stress_period_data={0:ssm_cells,1:ssm_cells,2:ssm_cells})

    nstrm = np.abs(m.sfr.nstrm)
    flopy.mt3d.Mt3dSft(mt,nsfinit=nstrm,mxsfbc=0,ietsfr=0,ioutobs=1001,nobssf=nstrm,
                       obs_sf=np.arange(nstrm)+1)

    mt.write_input()
    mt.run_model()




def setup_pest():
    m = flopy.modflow.Modflow.load(mf_nam, model_ws="template", check=False, exe_name=mf_exe)

    props = [["upw.hk",0],["rch.rech",0]]
    kperk = [[m.nper-1,0]]
    ph = pyemu.helpers.PstFromFlopyModel(mf_nam,org_model_ws="temp",new_model_ws="template",grid_props=props,
                                         const_props=props,sfr_pars=True,all_wells=True,remove_existing=True,
                                         model_exe_name=mf_exe,hds_kperk=kperk,
                                         extra_post_cmds=["{0} {1}".format(mt_exe,mt_nam)])


    m.external_path = '.'
    m.run_model()

    mt = flopy.mt3d.Mt3dms.load(mt_nam, model_ws="temp", exe_name=mt_exe)
    mt.change_model_ws("template")
    #mt.external_path = '.'
    mt.write_input()
    mt.run_model()

    df = write_ssm_tpl()

    os.chdir("template")

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

    ph.pst.parameter_groups.loc["kg","inctyp"] = "absolute"
    ph.pst.parameter_groups.loc["kg","derinc"] = derinc

    obs = ph.pst.observation_data
    obs.loc[:,"weight"] = 0.0
    sw_conc_obs = obs.loc[obs.obgnme=="sfrc","obsnme"]

    obs.loc[sw_conc_obs,"obgnme"] = "less_swconc"
    #obs.loc[sw_conc_obs, "weight"] = 1.0
    #only turn on one constraint in the middle of the domain
    obs.loc["sfrc20_1_03650.00","weight"] = 1.0

    # fix all non dec vars for now
    par = ph.pst.parameter_data
    par.loc[par.pargp!="kg","partrans"] = "fixed"

    ph.pst.pestpp_options = {}
    ph.pst.pestpp_options["opt_dec_var_groups"] = ["kg"]
    ph.pst.pestpp_options["opt_obj_func"] = 'obj.dat'
    ph.pst.pestpp_options["opt_direction"] = "max"
    ph.pst.pestpp_options["opt_risk"] = 0.5
    ph.pst.write(os.path.join("template","freyberg.pst"))
    pyemu.helpers.run("pestpp freyberg.pst",cwd="template")


    with open(os.path.join("template","obj.dat"),'w') as f:
        for pname in par.loc[par.pargp=="kg","parnme"]:
            f.write("{0} {1}\n".format(pname,1.0))

def write_ssm_tpl():
    ssm = os.path.join("template","freyberg.mt3d.ssm")
    f_in = open(ssm,'r')
    f_tpl = open(ssm+".tpl",'w')
    f_tpl.write("ptf ~\n")
    parnme = []
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
                pn = "~k_{0:02d}_{1:02d}~".format(r-1,c-1)
                line = " {0:9d} {1:9d} {2:9d} {3:9s} {4:9d}\n".format(l,r,c,pn,15)
                f_tpl.write(line)
                parnme.append(pn)
            f_tpl.write("-1\n-1\n-1\n-1\n-1\n")
            break
    df = pd.DataFrame({"tpl_str":parnme})
    df.loc[:,'parnme'] = df.tpl_str.apply(lambda x: x.replace("~",""))
    df.index = df.parnme
    df.loc[:,"pargp"] = "kg"
    df.loc[:,"parubnd"] = 10.0
    df.loc[:,"partrans"] = "none"
    df.loc[:,"parlbnd"] = 0.0
    df.loc[:,"parval1"] = 0.0
    return df


def run_jco():
    pst = pyemu.Pst(os.path.join("template","freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join("template","freyberg.pst"))

    pyemu.helpers.start_slaves("template","pestpp","freyberg.pst",num_slaves=15,slave_root='.',
                               master_dir="resp_master")

def run_pestpp_opt():
    pst_file = os.path.join("template", "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.control_data.noptmax = 3
    pst.write(pst_file)
    pyemu.helpers.start_slaves("template","pestpp-opt","freyberg.pst",num_slaves=15,master_dir="master_opt",
                               slave_root='.')

def spike_test():
    pst_file = os.path.join("template", "freyberg.pst")
    pst = pyemu.Pst(pst_file)
    pst.parameter_data.loc["k_10_00","parval1"] = derinc
    pst.write(pst_file)
    pyemu.helpers.run("pestpp {0}".format(os.path.split(pst_file)[-1]),cwd="template")
    pst = pyemu.Pst(pst_file)
    pst.res.loc[:,"wr"] = pst.res.residual * pst.res.weight
    print(pst.res.loc[pst.res.wr>0,:])
    print(pst.res.loc[pst.nnz_obs_names,:])

    m = flopy.modflow.Modflow.load(mf_nam, model_ws="template", check=False, exe_name=mf_exe)

    obs = pst.observation_data
    ucn_obs = obs.loc[obs.obsnme.apply(lambda x: x.startswith("ucn")),:]
    ucn_obs.loc[:,"i"] = ucn_obs.obsnme.apply(lambda x: int(x.split("_")[2]))
    ucn_obs.loc[:, "j"] = ucn_obs.obsnme.apply(lambda x: int(x.split("_")[3]))
    arr = np.zeros((m.nrow,m.ncol))
    arr[ucn_obs.i,ucn_obs.j] = pst.res.loc[ucn_obs.obsnme,"modelled"]
    arr = np.ma.masked_where(m.bas6.ibound[0].array==0,arr)
    plt.imshow(arr)
    plt.show()


if __name__ == "__main__":
    #setup_models()
    setup_pest()
    #write_ssm_tpl()
    #run_test()
    #spike_test()
    run_pestpp_opt()