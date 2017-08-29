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
    risk_vals.extend([x for x in np.arange(0.05,1.0,0.05)])
    risk_vals.extend([0.99])
    #risk_vals = [0.1,0.5,0.9]
    print(risk_vals)

    pst = pyemu.Pst("supply2_pest.fosm.pst")
    #pst.control_data.noptmax = 2
    # jco = pyemu.Jco.from_binary("supply2_pest.full.jcb")
    # par = pst.parameter_data
    # fosm_names = list(par.loc[par.pargp.apply(lambda x: x in ["trans","sfr_seg"]),"parnme"].values)
    # jco = jco.get(col_names=fosm_names)
    # jco.to_binary("supply2_pest.fosm.jcb")
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
        shutil.copy2("risk_sweep.1.par", os.path.join(results_dir, fname + ".par"))
        #break
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
        shutil.copy2("risk_sweep.1.par", os.path.join(results_dir, fname + ".par"))
        
    os.chdir("..")    

def check_infeas(f):
    with open(f,'r') as f:
        for line in f:
            if "infeasibility report" in line:
                return True
    return False

def plot_tradeoff():
   
   
    fig = plt.figure(figsize=(19.0 / 2.54, 9.0 / 2.54))
    ax = plt.subplot(111)

    rdir = os.path.join(wdir_base,results_dir)
    rei_files = [f for f in os.listdir(rdir) if f.endswith(".rei")]
    print(rei_files)
    dfs = [pyemu.pst_utils.read_resfile(os.path.join(rdir,f)) for f in rei_files]
    risk_vals = np.array([float(f.split('_')[1]) for f in rei_files])
    obj_func = np.array([df.loc["pi_obj_func","modelled"] for df in dfs])
    infeas = np.array([check_infeas(os.path.join(rdir,f.replace(".rei",".rec"))) for f in rei_files])
    obj_func[infeas==True] = np.NaN 
    ax.plot(risk_vals,obj_func,'b',lw=0.5,marker=".")


    # rdir = os.path.join(wdir_worth,results_dir)
    # rei_files = [f for f in os.listdir(rdir) if f.endswith(".rei")]
    # print(rei_files)
    # dfs = [pyemu.pst_utils.read_resfile(os.path.join(rdir,f)) for f in rei_files]
    # risk_vals = np.array([float(f.split('_')[1]) for f in rei_files])
    # obj_func = np.array([df.loc["pi_obj_func","modelled"] for df in dfs])
    # infeas = np.array([check_infeas(os.path.join(rdir,f.replace(".rei",".rec"))) for f in rei_files])
    # obj_func[infeas==True] = np.NaN 
    # ax.plot(risk_vals,obj_func,'b',lw=0.5,marker=".")

    # risk_infeas = risk_vals.copy()
    # risk_infeas[infeas==True] = np.NaN
    # first_infeas = np.nanmax(risk_infeas)
    # xlim = (0.0,1.0)
    ylim = ax.get_ylim()
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

    ax.plot((0.5,0.5),(ylim),color='r',lw=1.5)
    t = ax.text(0.5,ylim[1]*0.99,"risk neutral",ha='center',va="top",fontsize=10,rotation=90.0,color='r')
    t.set_bbox({"color":"w"})
    ypos = ylim[1]*0.98
    plt.annotate(s='', xy=(0.0,ypos), xytext=(0.475,ypos), arrowprops=dict(arrowstyle='<->',color='r'))
    t = ax.text(0.25,ypos,"risk tolerant",ha='center',va="center",fontsize=10,color='r')
    t.set_bbox({"color":"w"})
    plt.annotate(s='', xy=(0.525,ypos), xytext=(1.0,ypos), arrowprops=dict(arrowstyle='<->',color='r'))
    
    t = ax.text(0.75,ypos,"risk averse",ha='center',va="center",fontsize=10,color='r')
    t.set_bbox({"color":"w"})

    ax.set_xlabel("risk")
    ax.set_ylabel("optimal objective function value ($)")
    plt.tight_layout()
    #plt.show()
    plt.savefig("risk_tradeoff.pdf")

def get_gwm_decvar_vals(filename):
    dv_vals = {}
    with open(filename,'r') as f:
        while True:
            line = f.readline()
            if line == '':
                raise Exception()
            if "OPTIMAL RATES FOR EACH FLOW VARIABLE" in line:
                [f.readline() for _ in range(4)]
                while "----" not in line:
                    line = f.readline()
                    if"----" in line:
                        break
                    #print(line)
                    raw = line.lower().strip().split()
                    dv_vals[raw[0]] = float(raw[1])
                [f.readline() for _ in range(7)]
                while True:
                    line = f.readline()
                    if "----" in line:
                        break
                    raw = line.lower().strip().split()
                    dv_vals[raw[0]] = float(raw[2])
                
                return dv_vals




def plot_dev_var_bar():

    #parse the gwm output file for dec var values
    gwm_dvs = get_gwm_decvar_vals(os.path.join("baseline_opt","supply2.gwmout"))
    fig = plt.figure(figsize=(190.0 / 25.4, 70.0 / 25.4))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)


    rdir = os.path.join(wdir_base,results_dir)
    pst = pyemu.Pst(os.path.join(base,"supply2_pest.fosm.pst"))
    dvs = pst.parameter_data.loc[pst.parameter_data.pargp.apply(lambda x: x in ["pumping","external"]),"parnme"]
    par_files = [f for f in os.listdir(rdir) if f.endswith(".par")]
    print(par_files)
    dfs = [pd.read_csv(os.path.join(rdir,f),skiprows=1,header=None,names=["parnme","parval1"],
        usecols=[0,1],delim_whitespace=True,index_col=0).loc[dvs] for f in par_files]
    print(dfs[0])
    risk_vals = np.array([float(f.split('_')[1]) for f in par_files])
    # obj_func = np.array([df.loc["pi_obj_func","modelled"] for df in dfs])
    infeas = np.array([check_infeas(os.path.join(rdir,f.replace(".par",".rec"))) for f in par_files])
    print(infeas)
    # obj_func[infeas==True] = np.NaN 
    # ax.plot(risk_vals,obj_func,'b',lw=0.5,marker=".")
    nu_idx = np.argwhere(risk_vals==0.5)
    in_idx = np.argwhere(infeas==1)[0] - 1
    
    labels = ["A) risk = {0}".format(risk_vals[0]), "B) risk = 0.5", "C) risk = {0}".format(risk_vals[in_idx][0])]

    for ax, df, l  in zip([ax1,ax2,ax3],[dfs[0],dfs[nu_idx],dfs[in_idx]],labels):
        if "0.5" in l:
            print(gwm_dvs)
            df.loc[:,"gwm"] = df.index.map(lambda x: gwm_dvs[x])
            df.plot(kind="bar",ax=ax,legend=False,alpha=0.5)
        else:
            df.plot(kind="bar",ax=ax,legend=False,alpha=0.5)
        ax.set_xlabel("decision variable")
        ax.text(0.01,1.01,l,transform=ax.transAxes)
        


    ax1.set_ylabel("pumping rate ($\\frac{m^3}{d}$)")
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    for ax in [ax1,ax2,ax3]:
        ax.set_ylim(0,55000)
        ax.grid()
    plt.tight_layout()
    #plt.show()
    plt.savefig("dec_vars.pdf")



if __name__ == "__main__":
    #run_base()
    #run_worth()
    #plot_tradeoff()
    plot_dev_var_bar()