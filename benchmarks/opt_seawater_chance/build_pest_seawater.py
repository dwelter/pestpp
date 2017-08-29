import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu

plt_dir = "plots"
if not os.path.exists(plt_dir):
    os.mkdir(plt_dir)

nam = "SEAWATER.NAM"
org_model_ws = "SEAWATER"
new_model_ws = "template"
nlay,nrow,ncol = 2,20,30
delt = 200
hds_file = "seawater.hds"
grad_file = "gradients.dat"

def setup_model():

    om = flopy.modflow.Modflow.load(nam,model_ws=org_model_ws,check=False)
    om.change_model_ws(new_model_ws, reset_external=True)
    om.name = om.name.lower()
    om.external_path = "."
    flopy.modflow.ModflowOc(om,chedfm="({0}E15.6)".format(om.ncol),unitnumber=[50, 51, 52, 53, 0],
                            stress_period_data={(0,0):["save head","save budget"]},label='')
    om.remove_external(om.external_fnames[0])
    om.rch.rech[0].format.free = True

    om.write_input()
    om.run_model()
    mfl = flopy.utils.MfListBudget(os.path.join(om.model_ws,om.name+".list"))
    df_flx,df_vol = mfl.get_dataframes(diff=True)
    df_flx.iloc[0,:].plot(kind='bar')
    plt.savefig(os.path.join(plt_dir,"org_flx.png"))
    df_vol.iloc[0,:].plot(kind='bar')
    plt.savefig(os.path.join(plt_dir,"org_vol.png"))


def setup_wel_tpl():
    ifmt = lambda x: "{0:10d}".format(x)
    ffmt = lambda x: "{0:15.6E}".format(x)
    afmt = lambda x: " #"+x
    df = pd.read_csv("wel.csv")
    fmt = {"l":ifmt,"r":ifmt,"c":ifmt,"flux":ffmt,"aux":afmt}


    wel_file = os.path.join(new_model_ws,nam.lower().replace(".nam",".wel"))
    wel_name = os.path.join(new_model_ws, "WEL_0000.dat")
    tpl_name = wel_name + ".tpl"

    f_w = open(wel_file,'w')
    #mxwel
    f_w.write('{0:10d}{1:10d}\n'.format(13,0))
    f_w.write('{0:10d}{1:10d} \n        '.format(13, 0))
    f_w.write("open/close {0}\n".format(os.path.split(wel_name)[-1]))
    f_w.close()

    f_in = open(wel_name,'w')
    f_tpl = open(tpl_name, 'w')
    f_tpl.write("ptf ~\n")

    f_in.write("        ")
    f_tpl.write("        ")

    df.to_string(f_in,header=False,index=False,formatters=fmt)
    df.loc[:,"tpl_str"] = df.aux.apply(lambda x: "~    {0}     ~".format(x.split('_')[0]))
    df.loc[:,["l","r","c","tpl_str","aux"]].\
        to_string(f_tpl,header=False,index=False,formatters=fmt)
    f_in.close()
    f_tpl.close()

    m = flopy.modflow.Modflow.load(nam.lower(),model_ws=new_model_ws)
    m.run_model()
    mfl = flopy.utils.MfListBudget(os.path.join(m.model_ws, m.name + ".list"))
    df_flx, df_vol = mfl.get_dataframes(diff=True)
    df_flx.iloc[0, :].plot(kind='bar')
    plt.savefig(os.path.join(plt_dir, "mod_flx.png"))
    df_vol.iloc[0, :].plot(kind='bar')
    plt.savefig(os.path.join(plt_dir, "mod_vol.png"))


def setup_par_tpls():
    bcf_file = os.path.join(new_model_ws,nam.lower().replace(".nam",".bcf"))
    rch_file = bcf_file.replace(".bc6",".rch")
    for f in [bcf_file,rch_file]:
        f_in = open(f,'r')
        f_tpl = open(f+".tpl",'w')
        f_tpl.write("ptf ~\n")
        for line in f_in:

            if "open/close" in line.lower() and "wetdry" not in line.lower() \
                    and "factor" not in line.lower():
                raw = line.strip().split()
                raw[2] =  " ~    {0}     ~ ".format(raw[-1])
                line = ' '.join(raw) + '\n'
            f_tpl.write(line)

def setup_par_tpls():
    arr_files = ["rech_0.ref","Horizontal_Hydraulic_Conductivity_Layer_1.ref",
                 "Transmissivity_Layer_2.ref","Vertical_Conductance_Layer_1.ref"]
    prefixes = ["rch","hy1","tr2","vcd"]
    for arr_file,prefix in zip(arr_files,prefixes):
        with open(os.path.join(new_model_ws,arr_file+".tpl"),'w') as f:
            f.write("ptf ~\n")
            for i in range(nrow):
                for j in range(ncol):
                    pname = "{0}_{1:02d}_{2:02d}".format(prefix,i,j)
                    f.write(" ~   {0}    ~ ".format(pname))
                f.write("\n")


def setup_hds_obs():
    os.chdir(new_model_ws)
    process_hds()
    os.chdir("..")
    grad_df = pd.read_csv(os.path.join(new_model_ws,"gradients.dat"),delim_whitespace=True)
    grad_df.index = grad_df.obsnme
    grad_df.loc[:,"ins_str"] = grad_df.obsnme.apply(lambda x: "l1 w !{0}!".format(x))
    with open(os.path.join(new_model_ws,grad_file+".ins"),'w') as f:
        f.write("pif ~\nl1\n")
        grad_df.loc[:,["ins_str"]].to_string(f,index=False,header=False)
        f.write("\n")

    with open(os.path.join(new_model_ws,hds_file+".ins"),'w') as f:
        f.write("pif ~\n")
        for k in range(nlay):
            for i in range(nrow):
                f.write("l1 ")
                for j in range(ncol):
                    oname = "wl_l{0}_r{1:02}_c{2:02d}".format(k+1,i+1,j+1)
                    f.write(" !{0}! ".format(oname))
                f.write("\n")

    #super hack alert
    f_in = open("build_pest_seawater.py",'r')
    f_out = open(os.path.join(new_model_ws,"forward_run.py"),'w')
    f_out.write('import os\nimport numpy as np\n')
    f_out.write("hds_file = '{}'\n".format(hds_file))
    f_out.write("nrow,ncol,nlay = {0},{1},{2}\n".format(nrow,ncol,nlay))
    f_out.write("delt = {0}\n".format(delt))
    f_out.write("grad_file = '{0}'\n".format(grad_file))
    while True:
        line = f_in.readline()
        if line == '':
            raise Exception()
        if line.startswith("def process_hds():"):
            f_out.write(line)
            line = f_in.readline()
            while not line.startswith("def"):
                f_out.write(line)
                line = f_in.readline()
            break
    f_out.write("if __name__ == '__main__':\n")
    f_out.write("    os.system('mf2005 {0}')\n".format(nam.lower()))
    f_out.write("    process_hds()\n")
    f_out.close()
    f_in.close()

def process_hds():

    arr = np.loadtxt(hds_file).reshape((nlay,nrow,ncol))
    hg1 = (arr[0,:,27] - arr[0,:,29]) / (2.0 * delt)
    hg2 = (arr[1,:,27] - arr[0,:,29]) / (2.0 * delt)
    vg1 = (arr[0,:,28] - arr[1,:,28])
    print(arr[0,:,28], arr[1,:,28])
    hg1_labels = ["hg1_r{0}".format(i+1) for i in range(nrow)]
    hg2_labels = ["hg2_r{0}".format(i+1) for i in range(nrow)]
    vg1_labels = ["vg1_r{0}".format(i+1) for i in range(nrow)]
    with open(grad_file,'w') as f:
        f.write("obsnme obsval\n")
        for arr,labels in zip([hg1,hg2,vg1],[hg1_labels,hg2_labels,vg1_labels]):
            for v,l in zip(arr,labels):
                f.write("{0} {1:15.6E}\n".format(l,v))


def setup_pest():
    wel_df = pd.read_csv("wel.csv")
    os.chdir(new_model_ws)
    tpl_files = [f for f in os.listdir('.') if f.endswith(".tpl")]
    in_files = [f.replace(".tpl",'') for f in tpl_files]
    ins_files = [f for f in os.listdir('.') if f.endswith(".ins")]
    out_files = [f.replace(".ins",'') for f in ins_files]

    pst = pyemu.Pst.from_io_files(tpl_files,in_files,ins_files,out_files)
    par = pst.parameter_data
    p = par.parnme.apply(lambda x: x.startswith("w"))
    par.loc[p,"parval1"] = 0.0
    par.loc[p, "partrans"] = "none"
    par.loc[p, "parubnd"] = 10000.
    par.loc[p, "parlbnd"] = 0.
    par.loc[p, "scale"] = -1.0
    par.loc[p,"parchglim"] = "relative"
    par.loc[p, "pargp"] = "extrac_dv"

    p = par.parnme.apply(lambda x: x.startswith("i"))
    par.loc[p, "parval1"] = 0.0
    par.loc[p, "partrans"] = "none"
    par.loc[p, "parubnd"] = 10000.
    par.loc[p, "parlbnd"] = 0.
    par.loc[p, "scale"] = 1.0
    par.loc[p, "parchglim"] = "relative"
    par.loc[p, "pargp"] = "inject_dv"

    par.loc["ex", "parval1"] = 5000.0
    par.loc["ex", "partrans"] = "none"
    par.loc["ex", "parubnd"] = 10000.
    par.loc["ex", "parlbnd"] = 0.
    par.loc["ex", "scale"] = -1.0
    par.loc["ex", "parchglim"] = "relative"
    par.loc["ex", "pargp"] = "existing"

    prefixes = ["rch", "hy1", "tr2", "vcd"]
    parvals = [0.002,5.0,800.0,0.05]
    for prefix,parval in zip(prefixes,parvals):
        p = par.parnme.apply(lambda x: x.startswith(prefix))
        par.loc[p, "parval1"] = parval
        par.loc[p, "partrans"] = "log"
        par.loc[p, "parubnd"] = parval * 2.0
        par.loc[p, "parlbnd"] = parval * 0.5
        par.loc[p, "scale"] = 1.0
        par.loc[p, "pargp"] = prefix

    pst._rectify_pgroups()
    pargp = pst.parameter_groups
    dv_pargp = pargp.pargpnme.apply(lambda x: x in ["extrac_dv", "inject_dv"])
    pargp.loc[dv_pargp, "inctyp"] = "absolute"
    pargp.loc[dv_pargp, "derinc"] = 500.0

    obs = pst.observation_data
    obs.loc[:,"weight"] = 0.0
    obs.loc[:,"obgnme"] = obs.obsnme.apply(lambda x: x.split('_')[0])

    #horizontal grad constraints
    hg_obs = obs.loc[obs.obsnme.apply(lambda x: x.startswith("hg")),:]
    hg_obs.loc[:,"row"] = hg_obs.obsnme.apply(lambda x: int(x.split('_')[-1][1:]))
    hg_odd_rows = hg_obs.loc[hg_obs.row.apply(lambda x: x%2 != 0),"obsnme"]
    obs.loc[hg_odd_rows,"weight"] = 1.0
    obs.loc[hg_odd_rows, "obgnme"] = "greater_hg"
    obs.loc[hg_odd_rows, "obsval"] = 0.00375

    # vertical grad constraints
    vg_obs = obs.loc[obs.obsnme.apply(lambda x: x.startswith("vg")), :]
    vg_obs.loc[:, "row"] = vg_obs.obsnme.apply(lambda x: int(x.split('_')[-1][1:]))
    vg_odd_rows = vg_obs.loc[vg_obs.row.apply(lambda x: x % 2 != 0), "obsnme"]
    obs.loc[vg_odd_rows, "weight"] = 1.0
    obs.loc[vg_odd_rows, "obgnme"] = "greater_vg"
    obs.loc[vg_odd_rows,"obsval"] = 0.05

    # water levels at injection wells
    i_wel_df = wel_df.loc[wel_df.aux.apply(lambda x: x.startswith('I')),:]
    i_wel_df = i_wel_df.loc[i_wel_df.l==1,:]
    i_wel_df.loc[:,"obsnme"] = i_wel_df.apply(lambda x: "wl_l{0}_r{1:02d}_c{2:02d}".format(x.l,x.r,x.c),axis=1)
    obs.loc[i_wel_df.obsnme,"weight"] = 1.0
    obs.loc[i_wel_df.obsnme,"obsval"] = 5.0
    obs.loc[i_wel_df.obsnme,"obgnme"] = "less_wl"

    ext_wels = list(par.loc[par.parnme.apply(lambda x: x.startswith('w')),"parnme"])
    inj_wels = list(par.loc[par.parnme.apply(lambda x: x.startswith('i')),"parnme"])
    dec_vars = ext_wels
    dec_vars.extend(inj_wels)
    pst.add_pi_equation(dec_vars,obs_group="greater_sum",pilbl="greater_sum",coef_dict={iw:-1.0 for iw in inj_wels})
    pst.add_pi_equation(dec_vars, obs_group="obj_func", pilbl="obj_func", coef_dict={iw: -1.0 for iw in inj_wels})

    par.loc[par.parnme.apply(lambda x: x not in dec_vars),"partrans"] = "fixed"

    pst.model_command = ["python forward_run.py"]

    pst.control_data.noptmax = 0
    pst.pestpp_options["opt_direction"] = "max"
    pst.pestpp_options["opt_constraint_groups"] = "greater_sum,greater_vg,less_wl,greater_hg"
    pst.pestpp_options["opt_obj_func"] = "obj_func"
    pst.pestpp_options["opt_coin_log"] = "true"
    pst.pestpp_options["opt_dec_var_groups"] = "extrac_dv,inject_dv"
    pst.write(os.path.join("seawater_pest_init.pst"))
    pst.control_data.noptmax = 4
    pst.write(os.path.join("seawater_pest.pst"))
    
    os.system("pestchek seawater_pest.pst")
    os.system("ipestpp seawater_pest_init.pst")
    os.chdir("..")


if __name__ == "__main__":
    setup_model()
    setup_wel_tpl()
    setup_par_tpls()
    setup_hds_obs()
    setup_pest()