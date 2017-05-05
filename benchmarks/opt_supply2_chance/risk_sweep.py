import os
import subprocess as sp
import numpy as np
import pandas as pd
import shutil
import matplotlib
font = {'size'   : 8}
matplotlib.rc('font', **font)
from matplotlib.patches import Rectangle as rect
import matplotlib.pyplot as plt
import pyemu


base = "template"
wdir_base = "risk_sweep_base"
wdir_worth = "risk_sweep_worth"



results_dir = "results"

def run_base():
    if os.path.exists(wdir_base):
        shutil.rmtree(wdir_base)

    shutil.copytree(base,wdir_base)
    os.chdir(wdir_base)


    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    os.mkdir(results_dir)

    [os.remove(f) for f in os.listdir('.') if f.endswith(".rei")] 

    risk_vals = [0.01]
    risk_vals.extend([x for x in np.arange(0.025,1.0,0.025)])
    risk_vals.extend([0.99])
    #risk_vals = [0.1,0.5,0.9]
    print(risk_vals)

    pst = pyemu.Pst("supply2_pest.fosm.pst")

    for risk_val in risk_vals:
        pst.pestpp_options["opt_risk"] = risk_val
        pst.write("risk_sweep.pst")
        sp.call(["pestpp-opt.exe","risk_sweep.pst"])
        rei_files = [f for f in os.listdir('.') if f.endswith(".rei")]
        inums = [int(f.split('.')[-2]) for f in rei_files]
        mx_idx = inums.index(max(inums))
        fname = "risk_{0:04.3f}_".format(risk_val)
        shutil.copy2(rei_files[mx_idx],os.path.join(results_dir,fname+".rei"))
        shutil.copy2("risk_sweep.rec", os.path.join(results_dir, fname + ".rec"))
    os.chdir("..")


def run_worth():
    if os.path.exists(wdir_worth):
        shutil.rmtree(wdir_worth)

    shutil.copytree(base,wdir_worth)
    os.chdir(wdir_worth)

    worth_locs = [(10,20),(18,5),(8,15),(14,14)]

    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    os.mkdir(results_dir)

    [os.remove(f) for f in os.listdir('.') if f.endswith(".rei")] 

    risk_vals = [0.01]
    risk_vals.extend([x for x in np.arange(0.025,1.0,0.025)])
    risk_vals.extend([0.99])
    #risk_vals = [0.1,0.5,0.9]
    print(risk_vals)

    pst = pyemu.Pst("supply2_pest.fosm.pst")

    obs = pst.observation_data
    wl_obs = obs.loc[obs.obgnme=="head",:]
    wl_obs.loc[:,"i"] = wl_obs.obsnme.apply(lambda x: int(x[1:3]))
    wl_obs.loc[:,"j"] = wl_obs.obsnme.apply(lambda x: int(x.split('_')[-1]))
    wl_obs.loc[:,"ij"] = wl_obs.apply(lambda x: (x.i,x.j),axis=1)
    worth_names = wl_obs.loc[wl_obs.ij.apply(lambda x: x in worth_locs),"obsnme"]
    print(worth_names)
    pst.observation_data.loc[worth_names,"weight"] = 1.0

    for risk_val in risk_vals:
        pst.pestpp_options["opt_risk"] = risk_val
        pst.write("risk_sweep.pst")
        sp.call(["pestpp-opt.exe","risk_sweep.pst"])
        rei_files = [f for f in os.listdir('.') if f.endswith(".rei")]
        inums = [int(f.split('.')[-2]) for f in rei_files]
        mx_idx = inums.index(max(inums))
        fname = "risk_{0:04.3f}_".format(risk_val)
        shutil.copy2(rei_files[mx_idx],os.path.join(results_dir,fname+".rei"))
        shutil.copy2("risk_sweep.rec", os.path.join(results_dir, fname + ".rec"))
    os.chdir("..")    

def plot():

    def check_infeas(f):
        with open(os.path.join(rdir,f),'r') as f:
            for line in f:
                if "infeasibility report" in line:
                    return True
        return False
   
    fig = plt.figure(figsize=(19.0 / 2.54, 9.0 / 2.54))
    ax = plt.subplot(111)

    rdir = os.path.join(wdir_base,results_dir)
    rei_files = [f for f in os.listdir(rdir) if f.endswith(".rei")]
    print(rei_files)
    dfs = [pyemu.pst_utils.read_resfile(os.path.join(rdir,f)) for f in rei_files]
    risk_vals = np.array([float(f.split('_')[1]) for f in rei_files])
    obj_func = np.array([df.loc["pi_obj_func","modelled"] for df in dfs])
    infeas = np.array([check_infeas(f.replace(".rei",".rec")) for f in rei_files])
    obj_func[infeas==True] = np.NaN 
    ax.plot(risk_vals,obj_func,'b',lw=0.5,marker=".")


    rdir = os.path.join(wdir_worth,results_dir)
    rei_files = [f for f in os.listdir(rdir) if f.endswith(".rei")]
    print(rei_files)
    dfs = [pyemu.pst_utils.read_resfile(os.path.join(rdir,f)) for f in rei_files]
    risk_vals = np.array([float(f.split('_')[1]) for f in rei_files])
    obj_func = np.array([df.loc["pi_obj_func","modelled"] for df in dfs])
    infeas = np.array([check_infeas(f.replace(".rei",".rec")) for f in rei_files])
    obj_func[infeas==True] = np.NaN 
    ax.plot(risk_vals,obj_func,'b',lw=0.5,marker=".")

    # risk_infeas = risk_vals.copy()
    # risk_infeas[infeas==True] = np.NaN
    # first_infeas = np.nanmax(risk_infeas)
    # xlim = (0.0,1.0)
    # ylim = ax.get_ylim()
    # inf_rect = rect((first_infeas,ylim[0]),xlim[1]-first_infeas,ylim[1]-ylim[0],facecolor="c",alpha=0.5,edgecolor='none')
    # print(risk_vals)
    # ax.add_patch(inf_rect)
    xlim = (0.0,1.0)
    ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0.0,1.1,0.1))
    ax.grid()
    #t = ax.text((first_infeas + 1.0)/2.0,sum(ylim)/2.0,"infeasible region",ha='center',va="center",fontsize=10)
    #t.set_bbox({"color":"w"})
    ax.set_xlabel("risk")
    ax.set_ylabel("optimal objective function value ($)")
    plt.tight_layout()
    #plt.show()
    plt.savefig("risk_tradeoff.pdf")


if __name__ == "__main__":
    #run_base()
    #run_worth()
    plot()