import os
import string
import shutil

fontsize = 7

import matplotlib
font = {'size': fontsize}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import flopy
import pyemu

nrow,ncol = 40,20
num_reals = [30, 50, 100]
#num_reals = [30]
noptmax = 10
pst = pyemu.Pst(os.path.join("template", "pest.pst"))
pst.pestpp_options = {}
pst.pestpp_options["parcov"] = "prior.jcb"
pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]
pst.pestpp_options["lambda_scale_fac"] = 1.0
pst.pestpp_options["ies_subset_size"] = 10
pst.pestpp_options["ies_initial_lambda"] = 1.0
pst.observation_data.loc[pst.nnz_obs_names,"weight"] = 0.5

#pst.pestpp_options["parcov_filename"] = "prior.jcb"
#pst.parameter_data.loc[pst.adj_par_names,"partrans"] = "none"
#forecast_names = ["sw_gw_1","h01_28_11"]
forecast_names = ["c001fr35c11_19750101","flx_river_l_19750102","travel_time"]
ftitles = ["GW level","SW-GW exchange","travel time"]
cd = {30:'b',50:'g',100:'m'}
# pst.write(os.path.join("template","pest.pst"))
phi_accept = 13.0
m = flopy.modflow.Modflow.load("freyberg.nam", model_ws="template", check=False, forgive=False)
truth = np.loadtxt(os.path.join("template", "hk.truth.ref"))

full_fig = 190.0 / 25.4
onecol_fig = 90.0 / 25.4
onehalf_fig = 140.0 / 25.4



def run_pestpp():
    pst.control_data.noptmax = noptmax
    pst_name = "pest_base.pst"
    pst.write(os.path.join("template", pst_name))
    pyemu.helpers.start_slaves("template", "pestpp", pst_name,
                               num_slaves=10, master_dir="master_pestpp")

def run():
    for nr in num_reals:
    #for nr in [13]:
        pst.pestpp_options["ies_num_reals"] = nr
        # pst.pestpp_options["ies_use_prior_scaling"] = "true"
        # pst.pestpp_options["ies_initial_lambda"] = 1000000.0
        # master_dir = "master_{0}_ps".format(nr)
        # if os.path.exists(master_dir):
        #     shutil.rmtree(master_dir)
        # pst_name = "pest_{0}_ps.pst".format(nr)
        # pst.control_data.noptmax = noptmax
        # pst.write(os.path.join("template", pst_name))
        # pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
        #                            num_slaves=10, master_dir=master_dir)

        # pst.pestpp_options["ies_use_prior_scaling"] = "false"
        # pst.pestpp_options["ies_group_draws"] = "false"
        # pst.pestpp_options["ies_initial_lambda"] = 1.0
        # pst.pestpp_options["par_sigma_range"] = 100.0
        pst.control_data.noptmax = noptmax
        master_dir = "master_{0}_nps".format(nr)
        if os.path.exists(master_dir):
            shutil.rmtree(master_dir)
        pst_name = "pest_{0}_nps.pst".format(nr)
        pst.write(os.path.join("template", pst_name))
        pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
                                   num_slaves=10, master_dir=master_dir)


def plot_hk(hk_arr=None, ax=None,vmin=None,vmax=None):
    mm = flopy.plot.ModelMap(m.sr, ax=ax, model=m)
    save = True
    if hk_arr is None:
        save = False
        hk_arr = m.lpf.hk[0].array
    if vmin is None:
        vmin = hk_arr.min()
    if vmax is None:
        vmax = hk_arr.max()
    c = mm.plot_array(hk_arr, alpha=0.5, vmin=vmin, vmax=vmax)
    #if save:
    #    a = plt.colorbar(c)
    #    a.set_label("hk ($\\frac{m}{d}$)")
    mm.plot_ibound(m.bas6.ibound[0].array)
    mm.plot_bc("stage", package=m.riv)
    mm.plot_bc("flux", package=m.wel)

    # pi,pj = 14,2
    #pi = int(forecast_names[1].split('_')[1])
    #pj = int(forecast_names[1].split('_')[2])

    # head forecast
    pi,pj = 34,10
    px, py = m.sr.xcentergrid[pi, pj], m.sr.ycentergrid[pi, pj]
    ax.scatter([px], [py], marker='^', s=30, color='k')

    #particle starting location
    pi, pj = 14, 2
    px, py = m.sr.xcentergrid[pi, pj], m.sr.ycentergrid[pi, pj]
    ax.scatter([px], [py], marker='+', s=30, color='k')

    df = pst.observation_data.loc[pst.observation_data.obgnme=="calhead",:].copy()
    #df.loc[:, "x"] = df.apply(lambda x: m.sr.xcentergrid[x.row - 1, x.col - 1], axis=1)
    #df.loc[:, "y"] = df.apply(lambda x: m.sr.ycentergrid[x.row - 1, x.col - 1], axis=1)
    df.loc[:,'i'] = df.obsnme.apply(lambda x: int(x[6:8]))
    df.loc[:, 'j'] = df.obsnme.apply(lambda x: int(x[9:11]))
    #df.loc[:, "x"] = list(df.apply(lambda x: m.sr.xcentergrid[x.i, x.j], axis=1))
    y = list(df.apply(lambda x: m.sr.ycentergrid[x.i, x.j], axis=1))
    x = list(df.apply(lambda x: m.sr.xcentergrid[x.i, x.j], axis=1))

    ax.scatter(x, y, marker='.', color='k', s=20)
    # m.lpf.hk[0].plot(axes=[ax],colorbar=True)

    # rarr = m.riv.stress_period_data.array["stage"]
    # ax.imshow(m.sr.xcentergrid,m.sr.ycentergrid,rarr,color="r")

    # print(rarr)
    ax.set_xlabel("x (m)",labelpad=0.1)
    ax.set_ylabel("y (m)",labelpad=0.1)
    return c

def plot_domain():
    fig = plt.figure(figsize=(onecol_fig,onecol_fig * 1.25))
    ax = plt.subplot(111, aspect="equal")
    c = plot_hk(truth, ax)
    cb = plt.colorbar(c)
    cb.set_label("K ($\\frac{m}{d}$)")
    plt.tight_layout()
    plt.savefig("freyberg_domain.pdf")
    plt.close(fig)


def plot_phi():
    #df_mc = pd.read_csv(os.path.join("master_mc","pest_mc.phi.actual.csv"))
    master_dirs = [d for d in os.listdir('.') if "master_" in d and "pestpp" not in d and "mc" not in d]
    fig = plt.figure(figsize=(onehalf_fig,onehalf_fig))
    gs = gridspec.GridSpec(3,3)
    ax = plt.subplot(gs[:-1,:])
    ax1 = plt.subplot(gs[-1,0])
    ax2 = plt.subplot(gs[-1,1])
    ax3 = plt.subplot(gs[-1,-1])
    # ax = plt.subplot2grid((3,3),(0,0),rowspan=3,colspan=2)
    # ax2 = plt.subplot2grid((3,3),(0,0))
    # ax3 = plt.subplot2grid((3, 3), (0, 1))
    # ax3 = plt.subplot2grid((3, 3), (0, 2))
    #ax = plt.subplot(111)
    #ax_ps = plt.subplot(132)
    #ax_nps = plt.subplot(133)
    #dfs = []

    ad = {30:ax1,50:ax2,100:ax3}
    bins = np.arange(0,100,10)
    #for master_dir in master_dirs:
    for num_real in num_reals:
        for master_dir in master_dirs:

            phi_file = [f for f in os.listdir(master_dir) if "phi.actual" in f][0]
            nr = int(master_dir.split('_')[1])
            ps = master_dir.split('_')[-1]

            if nr != num_real:
                continue
            if ps != "nps":
                continue
            df = pd.read_csv(os.path.join(master_dir,phi_file))

            ls = '-'
            label = "{0} reals, prior scaling".format(nr)
            #print(ps)
            if ps == "nps":
                ls = '--'
                label = "{0} reals".format(nr)
            c = cd[nr]
            mean = df.loc[:,"mean"]
            std = df.loc[:,"standard_deviation"]
            #ax.plot(df.total_runs,mean,ls=ls,color=c,lw=1.0,marker='.',ms=4,label=label)
            ax.semilogy(df.total_runs, mean, ls=ls, color=c, lw=1.0, marker='.', ms=4, label=label)
            #ax.fill_between(df.total_runs,mean-std,mean+std,facecolor=c,
            #                alpha=0.1,hatch='/////')
            #[ax.plot(df.total_runs,df.iloc[:,i],ls=ls,color=c,lw=0.05,label='') for i in range(6,df.shape[1])]
            #mask = df.iloc[-1,5:] > 500.0
            #df.iloc[-1,5:].loc[mask] = np.NaN
            if ps == "ps":
                df.iloc[-1,5:].hist(ax=ad[nr],facecolor=c,alpha=0.5,bins=bins,
                                    normed=True,edgecolor="none",)
            else:
                df.iloc[-1, 5:].hist(ax=ad[nr], bins=bins, normed=True,hatch='/////',
                                    edgecolor=c,facecolor="none")
            #df_mc.iloc[0,5:].hist(ax=ad[nr],bins=10,normed=True,facecolor='0.5',edgecolor="none")


    #ax.set_xlim(0,ax.get_xlim()[1])
    #ax.set_xlim(0,pst.npar_adj)
    #ax.set_xlim(0, 1000)

    #ax.set_ylim(0, 10000)
    ax.set_xlabel("model runs")
    ax.set_ylabel("objective function ($log_{10}\phi$)")
    ax.grid()
    ax.legend()
    ax.set_title('A)',loc="left",fontsize=fontsize)
    mx,mn = -1.0e+10,1.0e+10
    for a,t in zip([ax1,ax2,ax3],['B) 30 reals','C) 50 reals','D) 100 reals']):
        a.set_yticks([])
        a.set_xticks(bins[::2])
        a.set_xlabel('$\phi$')
        a.set_title(t,loc='left',fontsize=fontsize)
        #a.set_xlim(0,50)
        a.grid()
        xlim = a.get_xlim()
        mx = max(mx,xlim[1])
        mn = min(mn,xlim[0])
    for a in [ax1, ax2, ax3]:
        a.set_xlim(mn,mx)

    plt.tight_layout()
    plt.savefig("freyberg_phi.pdf")


def run_mc():
    pst.pestpp_options["ies_num_reals"] = 100000
    pst.control_data.noptmax = -1

    pst_name = "pest_mc.pst"
    pst.write(os.path.join("template",pst_name))
    master_dir = "master_mc"
    pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
                               num_slaves=10, master_dir=master_dir)


def plot_histograms():
    master_dirs = [d for d in os.listdir('.') if "master_" in d and "pestpp" not in d and "mc" not in d]
    df_mc_phi = pd.read_csv(os.path.join("master_mc", "pest_mc.phi.actual.csv")).iloc[:,6:].transpose()
    #print(df_mc_phi)
    keep_reals = [int(k) for k in list(df_mc_phi[(df_mc_phi.values < 21)].index)]
    print(len(keep_reals))
    df_mc = pd.read_csv(os.path.join("master_mc","pest_mc.0.obs.csv"))
    df_mc.columns = df_mc.columns.map(str.lower)
    #for f in forecast_names:
    #    print(df_mc.loc[keep_reals,f])
    #return
    fig = plt.figure(figsize=(full_fig,onehalf_fig))
    gs = gridspec.GridSpec(len(forecast_names),3)
    noptstr = ".{0}.".format(noptmax)
    lb = ["A)","B)","C)","D)","E)","F)"]
    axes = {f:[] for f in forecast_names}
    abet = string.ascii_uppercase
    for j,num_real in enumerate([30,50,100]):
        for master_dir in master_dirs:
            nr = int(master_dir.split('_')[1])
            ps = master_dir.split('_')[-1]
            if ps != "nps":
                continue
            if nr != num_real:
                continue
            print(master_dir,noptstr)
            ocsv_pt = [f for f in os.listdir(master_dir) if "obs.csv" in f and not "base" in f][-1]
            print(ocsv_pt)
            df_pt = pd.read_csv(os.path.join(master_dir,ocsv_pt))
            df_pt.columns = df_pt.columns.map(str.lower)
            for i,f in enumerate(forecast_names):
                ax = plt.subplot(gs[i,j])
                bins = np.linspace(df_mc.loc[keep_reals,f].min(),df_mc.loc[keep_reals,f].max(),10)
                df_mc.loc[keep_reals,f].hist(ax=ax,bins=bins,facecolor="0.5",edgecolor="none",alpha=0.5,normed=True)
                if ps == "ps":
                    df_pt.loc[:,f].hist(ax=ax,bins=bins,facecolor=cd[num_real],alpha=0.25,normed=True)
                else:
                    df_pt.loc[:, f].hist(ax=ax, bins=bins, facecolor='none',alpha=0.25,edgecolor=cd[num_real],hatch='////',normed=True)

                ax.set_yticks([])
                ax.grid()
                #ax.set_title(lb[j + (i*3)]+" {0}, {1} reals".format(f,nr),loc="left")
                ax.set_title(abet[j+(i+(len(forecast_names)*j))] + ") {0},{1} reals".\
                             format(ftitles[i], nr), loc="left",fontsize=fontsize)
                if f.startswith('travel'):
                    ax.set_xlabel("travel time ($d$)")
                elif f.startswith("c"):
                    ax.set_xlabel("gw level (m)",labelpad=0.1)
                else:
                    ax.set_xlabel("sw-gw exchange ($\\frac{m^3}{d}$)",labelpad=0.1)
                #ylim = ax.get_ylim()
                #v = pst.observation_data.loc[f,"obsval"]
                #ax.plot([v,v],ylim,'k',lw=2.0)
                #ax.set_ylim(ylim)
                axes[f].append(ax)

    for f, axs in axes.items():
        mn,mx = 1.0e+10,-1.0e10
        for ax in axs:
            n,x = ax.get_ylim()
            mn = min(mn,n)
            mx = max(mx,x)
        for ax in axs:
            ax.set_ylim(mn,mx)
            v = pst.observation_data.loc[f, "obsval"]
            ax.plot([v, v], [mn,mx], 'k', lw=2.0)

        # for f,ax in zip(forecast_names,axs):
        #     mc = df_mc.loc[keep_reals,f]
        #     print(mc)
        #     df_mc.loc[keep_reals,f].hist(bins=10,ax=ax,facecolor="0.5",normed=True)
        #     ax.set_xlim(mn,mx)
    plt.tight_layout()
    plt.savefig("freyberg_hist.pdf")
    plt.close(fig)


def get_hk_arr(df,real=None):
    if real is None:
        real = df.index[0]

    hk_pars = [c for c in df.columns.values if c.startswith('hk')]
    iidx = [int(x.split('_')[1][1:]) for x in hk_pars]
    jidx = [int(x.split('_')[2][1:]) for x in hk_pars]
    arr = np.zeros((nrow,ncol))
    arr[iidx,jidx] =df.loc[real,hk_pars]
    print(arr.shape)
    return arr



def plot_hk_arrays(real=None):
    fig = plt.figure(figsize=(full_fig,full_fig))
    gs = gridspec.GridSpec(2, 3,bottom=0.135)
    master_dirs = [d for d in os.listdir('.') if "master_" in d and "pestpp" not in d and "mc" not in d]
    noptstr = ".{0}.".format(noptmax)
    vmin = truth.min()
    vmax = truth.max()
    abet = string.ascii_uppercase
    c = 0
    for j,num_real in enumerate(num_reals):
        for master_dir in master_dirs:
            nr = int(master_dir.split('_')[1])
            ps = master_dir.split('_')[-1]
            if ps != "nps":
                continue
            if nr != num_real:
                continue

            pcsv_pt = [f for f in os.listdir(master_dir) if ".par" in f and noptstr in f][0]
            pcsv_pr = pcsv_pt.replace(noptstr,".0.")

            df_phi = pd.read_csv(os.path.join(master_dir,pcsv_pr.replace(".0.par",".phi.actual")))
            df_pt = pd.read_csv(os.path.join(master_dir,pcsv_pt),index_col=0)
            df_pr = pd.read_csv(os.path.join(master_dir, pcsv_pr),index_col=0)
            df_pt.columns = df_pt.columns.map(str.lower)
            df_pr.columns = df_pr.columns.map(str.lower)
            print(df_pr.columns )
            df_pr.loc["base",pst.par_names] = pst.parameter_data.parval1
            df_pr.index = [str(i) for i in df_pr.index]
            #print(df_pr.iloc[-1,:])
            #real_pr = list(df_pr.index).index("base")
            #real_pt = list(df_pt.index).index("base")

            if real is None:
                real = df_pt.index[0]
            real_phi_pt = df_phi.iloc[-1, :][real]
            real_phi_pr = df_phi.iloc[0, :][real]

            arr = get_hk_arr(df_pr,real)
            ax = plt.subplot(gs[0,j],aspect="equal")
            plot_hk(arr,ax,vmin=vmin,vmax=vmax)
            if j != 0:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.set_title("{0}) prior, {1} reals\nphi:{2}".format(abet[c],num_real,real_phi_pr),
                         fontsize=fontsize)
            c +=1

            arr = get_hk_arr(df_pt,real)
            ax = plt.subplot(gs[1, j],aspect="equal")
            cb = plot_hk(arr, ax, vmin=vmin, vmax=vmax)
            if j != 0:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            ax.set_title("{0}) post, {1} reals\nphi:{2}".format(abet[c], num_real,real_phi_pt),
                         fontsize=fontsize)
            c+= 1
            #break
    #fig = plt.figure(figsize=(6.5, 6.5))
    ax = plt.axes((0.1,0.06,0.8,0.02))
    cb = plt.colorbar(cb,cax=ax,orientation="horizontal")
    cb.set_label("K ($\\frac{m}{d}$)",labelpad=0.1)
    plt.tight_layout()
    plt.savefig("freyberg_par_{0}.pdf".format(real))
    plt.close(fig)

def plot_hk_arrays_figure():
    fig = plt.figure(figsize=(full_fig,full_fig))
    gs = gridspec.GridSpec(3, 3,bottom=0.135)
    master_dirs = [d for d in os.listdir('.') if "master_" in d and "pestpp" not in d and "mc" not in d]
    noptstr = ".{0}.".format(noptmax)
    vmin = truth.min()
    vmax = truth.max()
    abet = string.ascii_uppercase
    c = 0
    for j,num_real in enumerate(num_reals):
        for master_dir in master_dirs:
            nr = int(master_dir.split('_')[1])
            ps = master_dir.split('_')[-1]
            if ps != "nps":
                continue
            if nr != num_real:
                continue

            #pcsv_pt = [f for f in os.listdir(master_dir) if ".par" in f and noptstr in f][0]
            pcsv_pt = [f for f in os.listdir(master_dir) if ".par.csv" in f][-1]
            pcsv_pr = [f for f in os.listdir(master_dir) if ".par.csv" in f][0]

            df_phi = pd.read_csv(os.path.join(master_dir,pcsv_pr.replace(".0.par",".phi.actual")))
            print(df_phi.columns)
            df_pt = pd.read_csv(os.path.join(master_dir,pcsv_pt),index_col=0)
            df_pr = pd.read_csv(os.path.join(master_dir, pcsv_pr),index_col=0)
            df_pt.columns = df_pt.columns.map(str.lower)
            df_pr.columns = df_pr.columns.map(str.lower)
            print(df_pr.columns )
            df_pr.loc["base",pst.par_names] = pst.parameter_data.parval1
            df_pr.index = [str(i) for i in df_pr.index]
            #print(df_pr.iloc[-1,:])
            #real_pr = list(df_pr.index).index("base")
            #real_pt = list(df_pt.index).index("base")

            #if real is None:
            #    real = df_pt.index[0]
            real = '0'
            #print(df_phi.columns)
            real_phi_pt = df_phi.iloc[-1, :][real]
            real_phi_pr = df_phi.iloc[0, :][real]

            arr = get_hk_arr(df_pr,real)
            ax = plt.subplot(gs[0,j],aspect="equal")
            plot_hk(arr,ax,vmin=vmin,vmax=vmax)
            if j != 0:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.set_title("{0}) prior, {1} reals\nphi:{2:4G}".format(abet[c],num_real,real_phi_pr),
                         fontsize=fontsize)
            c +=1

            arr = get_hk_arr(df_pt,real)
            ax = plt.subplot(gs[1, j],aspect="equal")
            cb = plot_hk(arr, ax, vmin=vmin, vmax=vmax)
            if j != 0:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.set_title("{0}) posterior, {1} reals\nphi:{2:4G}".format(abet[c], num_real, real_phi_pt),
                         fontsize=fontsize)
            c+= 1

            real_phi_pt = df_phi.iloc[-1, :]["base"]
            arr = get_hk_arr(df_pt, "base")
            ax = plt.subplot(gs[2, j], aspect="equal")
            cb = plot_hk(arr, ax, vmin=vmin, vmax=vmax)
            if j != 0:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            ax.set_title("{0}) posterior, {1} reals\nphi:{2:4G}".format(abet[c], num_real, real_phi_pt),
                         fontsize=fontsize)
            c += 1

            #break
    #fig = plt.figure(figsize=(6.5, 6.5))
    ax = plt.axes((0.25,0.06,0.5,0.02))
    cb = plt.colorbar(cb,cax=ax,orientation="horizontal")
    cb.set_label("K ($\\frac{m}{d}$)",labelpad=0.1)
    plt.tight_layout()
    plt.savefig("freyberg_par_{0}.pdf".format(real))
    plt.close(fig)


def sigma_range_invest():
    nr = 13
    dfs,fore_dfs = {},{}
    sigma_ranges = [2,4,8,20]
    for sr in sigma_ranges:
        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = nr

        pst.pestpp_options["ies_use_prior_scaling"] = "false"
        pst.pestpp_options["ies_group_draws"] = "false"
        pst.pestpp_options["ies_initial_lambda"] = 1.0
        pst.pestpp_options["par_sigma_range"] = sr
        pst.control_data.noptmax = noptmax
        master_dir = "test"
        if os.path.exists(master_dir):
            shutil.rmtree(master_dir)
        pst_name = "pest_{0}_nps.pst".format(sr)
        pst.write(os.path.join("template", pst_name))
        pyemu.helpers.start_slaves("template", "pestpp-ies", pst_name,
                                   num_slaves=14, master_dir=master_dir)
        df = pd.read_csv(os.path.join(master_dir,pst_name.replace(".pst",".phi.actual.csv")),index_col=0)
        dfs[sr] = df
        fore_df = pd.read_csv(os.path.join(master_dir, pst_name.replace(".pst", ".{0}.obs.csv".format(noptmax))),
                              index_col=0)
        fore_df.columns = fore_df.columns.map(str.lower)
        fore_df = fore_df.loc[:, forecast_names]
        fore_dfs[sr] = fore_df
        #break


    fig = plt.figure(figsize=(20,20))
    ax = plt.subplot(111)
    colors = ['b','r','g','y','m','c','k','0.5']
    for i,(sr,df) in enumerate(dfs.items()):
        c = colors[i]

        for ii,rname in enumerate(df.columns[6:]):
            if ii == 0:
                label = sr
            else:
                label = None
            #ax.plot(df.total_runs,df.loc[:,rname].apply(np.log10).values,color=c,lw=0.5,label=label,alpha=0.5)
            ax.plot(df.total_runs,df.loc[:,"mean"].apply(np.log10).values,color=c,lw=0.5,label=label)
    ax.legend()
    plt.savefig("sigmarange.pdf")
    plt.close(fig)

    #df = pd.DataFrame(fore_dfs)

    std_dfs = {sr:df.std() for sr,df in fore_dfs.items()}
    std_df = pd.DataFrame(std_dfs)
    for fore in forecast_names:
        std_df.loc[fore,:].plot(kind="bar")
        plt.savefig("sigmarange_{0}.pdf".format(fore))
        plt.close("all")

    color_dfs = {colors[i]:fore_dfs[sr] for i,(sr) in enumerate(sigma_ranges)}
    pyemu.plot_utils.ensemble_helper(color_dfs,deter_vals=pst.observation_data.obsval.to_dict())

    # df.plot(kind="bar")
    # plt.savefig("sigmarange_forestd.pdf")
    # plt.close("all")
    print(df)

if __name__ == "__main__":
    #run()
    #sigma_range_invest()
    #run_pestpp()
    #run_mc()
    #plot_domain()
    #plot_phi()
    #plot_hk_arrays("base")
    #plot_hk_arrays()
    plot_histograms()
    #plot_hk_arrays_figure()
