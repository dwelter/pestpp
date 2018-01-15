import os
import shutil
import pandas as pd
import pyemu
pst_file = "pest.pst"
pst = pyemu.Pst(pst_file)


def gen_fails():
    mc = pyemu.MonteCarlo(pst=pst)
    mc.draw(5,obs=True)
    mc.parensemble.iloc[0,:] *= -10000000.0
    mc.parensemble.to_csv("par.csv")
    mc.obsensemble.to_csv("obs.csv")

    os.system("sweep.exe {0}".format(pst_file))
    df = pd.read_csv("sweep_out.csv")
    print(df.loc[:,pst.observation_data.obsnme.apply(str.upper)])


def gen_full():
    mc = pyemu.MonteCarlo(pst=pst)
    mc.draw(30,obs=True)
    mc.parensemble.to_csv("par.csv")
    mc.obsensemble.to_csv("obs.csv")
    mc.parensemble.to_csv("par1.csv")
    mc.obsensemble.to_csv("obs1.csv")

    os.system("sweep.exe {0}".format(pst_file))
    df = pd.read_csv("sweep_out.csv")
    df = df.loc[:,pst.observation_data.obsnme.apply(str.upper)]
    df.to_csv("obs_restart.csv")
    shutil.copy2("sweep_out.csv","sweep_out.bak.csv")

def ies():
    pst.control_data.noptmax = 1
    ies = pyemu.EnsembleSmoother(pst=pst,verbose="ies.log",save_mats=True)
    ies.initialize(parensemble="par1.csv",obsensemble="obs1.csv",restart_obsensemble="sweep_out.bak.csv")
    for i in range(6):
        ies.update(lambda_mults=[1.0],use_approx=False)

    #ies.update(lambda_mults=[10.,1,0.1],run_subset=10)
if __name__ == "__main__":
    gen_full()
    ies()