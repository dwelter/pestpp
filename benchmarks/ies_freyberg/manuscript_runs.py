import os
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import flopy
import pyemu

num_reals = [30, 50, 100]
noptmax = 10
pst = pyemu.Pst(os.path.join("template", "pest.pst"))


# pst.write(os.path.join("template","pest.pst"))
def run():
    for nr in num_reals:
        pst.pestpp_options["ies_num_reals"] = nr
        pst.pestpp_options["ies_use_prior_scaling"] = "true"
        master_dir = "master_{0}_ps".format(nr)
        if os.path.exists(master_dir):
            shutil.rmtree(master_dir)
        pst_name = "pest_{0}_ps.pst".format(nr)
        pst.control_data.noptmax = noptmax
        pst.write(os.path.join("template", pst_name))
        pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
                                   num_slaves=10, master_dir=master_dir)
        pst.pestpp_options["ies_use_prior_scaling"] = "false"
        master_dir = "master_{0}_nps".format(nr)
        if os.path.exists(master_dir):
            shutil.rmtree(master_dir)
        pst_name = "pest_{0}_nps.pst".format(nr)
        pst.write(os.path.join("template", pst_name))
        pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
                                   num_slaves=10, master_dir=master_dir)


def plot_domain():
    fig = plt.figure(figsize=(3.5, 4))
    ax = plt.subplot(111, aspect="equal")

    m = flopy.modflow.Modflow.load("freyberg.truth.nam", model_ws="template", check=False, forgive=False)
    mm = flopy.plot.ModelMap(m.sr,ax=ax,model=m)
    c = mm.plot_array(m.lpf.hk[0].array,alpha=0.5)
    a = plt.colorbar(c)
    a.set_label("hk ($\\frac{m}{d}$)")
    mm.plot_ibound(m.bas6.ibound[0].array)
    mm.plot_bc("stage",package=m.riv)
    mm.plot_bc("flux",package=m.wel)

    pi,pj = 14,2
    px,py = m.sr.xcentergrid[pi,pj],m.sr.ycentergrid[pi,pj]
    ax.scatter([px],[py],marker='+',s=30,color='k')

    df = pd.read_csv(os.path.join("template","obs_rowcol.dat"),delim_whitespace=True)
    df.loc[:,"x"] = df.apply(lambda x: m.sr.xcentergrid[x.row-1,x.col-1],axis=1)
    df.loc[:, "y"] = df.apply(lambda x: m.sr.ycentergrid[x.row - 1, x.col - 1], axis=1)

    ax.scatter(df.x,df.y,marker='.',color='k',s=20)
    #m.lpf.hk[0].plot(axes=[ax],colorbar=True)

    #rarr = m.riv.stress_period_data.array["stage"]
    #ax.imshow(m.sr.xcentergrid,m.sr.ycentergrid,rarr,color="r")

    #print(rarr)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plt.tight_layout()
    plt.savefig("freyberg_domain.pdf")
    plt.close(fig)


def plot_phi():
    master_dirs = [d for d in os.listdir('.') if "master_" in d and "pestpp" not in d]
    fig = plt.figure(figsize=(6.5,3.5))
    ax = plt.subplot(111)
    #ax_ps = plt.subplot(132)
    #ax_nps = plt.subplot(133)
    #dfs = []
    cd = {30:'b',50:'g',100:'m'}
    for master_dir in master_dirs:
        phi_file = [f for f in os.listdir(master_dir) if "phi.actual" in f][0]
        nr = int(master_dir.split('_')[1])
        ps = master_dir.split('_')[-1]
        df = pd.read_csv(os.path.join(master_dir,phi_file))
        ls = '-'
        label = "{0} reals, prior scaling".format(nr)
        print(ps)
        if ps == "nps":
            ls = '--'
            label = "{0} reals, no prior scaling".format(nr)


        c = cd[nr]
        ax.plot(df.total_runs,df.loc[:,"mean"],ls=ls,color=c,lw=1.0,marker='.',label=label)
        #[ax.plot(df.total_runs,df.iloc[:,i],ls=ls,color=c,lw=0.05) for i in range(5,df.shape[1])]
    ax.set_xlim(0,ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_xlabel("model runs")
    ax.set_ylabel("mean objective function ($\phi$)")
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.savefig("freyberg_phi.pdf")





def plot_histograms():
    pass


if __name__ == "__main__":
    # run()
    #plot_domain()
    plot_phi()
    plot_histograms()
