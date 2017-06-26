import os
import pandas as pd
import pyemu
pst_file = "pest.pst"
pst = pyemu.Pst(pst_file)


def gen():
    mc = pyemu.MonteCarlo(pst=pst)
    mc.draw(5,obs=True)
    mc.parensemble.iloc[0,:] *= -10000000.0
    mc.parensemble.to_csv("par.csv")
    mc.obsensemble.to_csv("obs.csv")

    os.system("sweep.exe {0}".format(pst_file))
    df = pd.read_csv("sweep_out.csv")
    print(df.loc[:,pst.observation_data.obsnme.apply(str.upper)])


def ies():
    ies = pyemu.EnsembleSmoother(pst=pst)
    ies.initialize(parensemble="par.csv",obsensemble="obs.csv")
    ies.update()

if __name__ == "__main__":
    #gen()
    ies()