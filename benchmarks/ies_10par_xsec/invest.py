import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pyemu
base_d = "template"
oname = "h01_06"

def setup():
    if os.path.exists("master"):
        shutil.rmtree("master")
    #shutil.copytree(base_d,"master")
    pst = pyemu.Pst(os.path.join(base_d,"pest.pst"))
    pst.pestpp_options = {}
    pst.control_data.noptmax = 3
    #pst.pestpp_options["ies_subset_size"] = 10

    nreal_per = 1

    dfs = []
    obsvals = np.linspace(0.0,5.0,20)

    for obsval in obsvals:
        obs = pst.observation_data.obsval.copy()
        obs.loc[oname] = obsval
        # for i in range(nreal_per):
        #     dfs.append(obs)
        for obsval2 in obsvals:
            obs1 = obs.copy()
            obs1.loc["h01_04"] = obsval2
            dfs.append(obs)

    df = pd.concat(dfs,axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(base_d,"obs.csv"))

    pst.pestpp_options["ies_obs_en"] = "obs.csv"
    pst.pestpp_options["ies_num_reals"] = df.shape[0]
    pst.write(os.path.join(base_d, "pest_pareto.pst"))

    #pyemu.os_utils.run("pestpp-ies pest.pst /h :4004",cwd="master")
    pyemu.os_utils.start_slaves("template","pestpp-ies","pest_pareto.pst",num_slaves=20,master_dir="master")

    plot()

def plot():
    pst = pyemu.Pst(os.path.join("master","pest_pareto.pst"))
    fig = plt.figure(figsize=(12, 12))
    ax = plt.subplot(221)
    ax2 = plt.subplot(224)
    ax3 = plt.subplot(223)

    df_init = df_init = pd.read_csv(os.path.join("master", "pest_pareto.0.obs.csv"), index_col=0)
    df_phi = pd.read_csv(os.path.join("master", "pest_pareto.phi.actual.csv")).loc[:, df_init.index]
    colors = ["0.5","b","c","b","r"]
    #for i in range(1,pst.control_data.noptmax):
    for i in [1,pst.control_data.noptmax]:
        df_final = pd.read_csv(os.path.join("master", "pest_pareto.{0}.obs.csv".format(pst.control_data.noptmax)), index_col=0)
        ax.scatter(df_phi.iloc[i-1,:], df_init.H01_06, color=colors[i-1],s=6)
        ax.scatter(df_phi.iloc[i,:], df_final.H01_06,color=colors[i],s=6)
        ax2.scatter(df_phi.iloc[i - 1, :], df_init.H01_04, color=colors[i - 1], s=6)
        ax2.scatter(df_phi.iloc[i, :], df_final.H01_04, color=colors[i], s=6)
        ax3.scatter(df_init.H01_06, df_init.H01_04, color=colors[i - 1], s=6)
        ax3.scatter(df_final.H01_06, df_final.H01_04, color=colors[i], s=6)


        for real in df_init.index:
            xi = df_phi.iloc[i-1,:][real]
            yi = df_init.loc[real,oname.upper()]
            xf = df_phi.iloc[i,:][real]
            yf = df_final.loc[real,oname.upper()]
            ax.plot([xi,xf],[yi,yf],color=colors[i-1],lw=0.15,zorder=0,alpha=0.5)

            xi = df_phi.iloc[i - 1, :][real]
            yi = df_init.loc[real, "H01_04"]
            xf = df_phi.iloc[i, :][real]
            yf = df_final.loc[real, "H01_04"]
            ax2.plot([xi, xf], [yi, yf], color=colors[i - 1], lw=0.15, zorder=0, alpha=0.5)

            xi = df_init.loc[real, "H01_06"]
            yi = df_init.loc[real, "H01_04"]
            xf = df_final.loc[real,"H01_06"]
            yf = df_final.loc[real, "H01_04"]
            ax3.plot([xi, xf], [yi, yf], color=colors[i - 1], lw=0.15, zorder=0, alpha=0.5)
        df_init = df_final
    ax.set_xlabel("phi")
    ax.set_ylabel(oname)
    ax2.set_xlabel("phi")
    ax2.set_ylabel("h01_04")

    ax3.set_xlabel("h01_06")
    ax3.set_ylabel("h01_04")

    ax.grid()
    ax2.grid()
    ax3.grid()

    plt.savefig("pareto.pdf")



if __name__ == "__main__":
    #setup()
    plot()