import os
import numpy as np
import pandas as pd

sfr_file = "supply2.sfrout"
aq_ex_name = "sfr_aq_ex.dat"
org_ex_name = "sfr_aq_ex.dat.org"
nrch = 40

def load_sfrout_aq_ex(filename):
    seg_idx,rech_idx = 3,4
    data_idx = 7 #stream flow
    #build a list of file starting positions
    starts = {}
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if line == '':
                break
            [f.readline() for _ in range(2)]
            line = f.readline().strip().split()
            kper,kstp = int(line[3]),int(line[5])
            [f.readline() for _ in range(4)]
            #kperkstp = "sp{0:02d}ts{1:02d}".format(kper+1,kstp+1)
            kperkstp = "{0:02d}".format(kper)
            starts[kperkstp] = f.tell()
            [f.readline() for _ in range(nrch)]

    dfs = []
    for kperkstp,pos in starts.items():
        with open(filename,'r') as f:
            f.seek(pos)
            df = pd.read_csv(f, delim_whitespace=True, nrows=nrch, usecols=[seg_idx, rech_idx, data_idx],
                             names=["segment", "reach", "aq_ex"])
        df.loc[:,kperkstp] = df.pop("aq_ex")
        df.loc[:,"segrech"] = df.segment.apply(lambda x: "s{0:01d}".format(x)) + \
                             df.reach.apply(lambda x: "r{0:02d}".format(x))
        df.pop("segment")
        df.pop("reach")
        df.index = df.pop("segrech")
        dfs.append(df)
    df = pd.concat(dfs,axis=1)
    cols = list(df.columns)
    cols.sort()
    return df.loc[:,cols]


def apply():
    df_org = pd.read_csv(org_ex_name,delim_whitespace=True,index_col=0)
    df = load_sfrout_aq_ex(sfr_file)
    delta_df = df_org - df
    delta_df.to_csv(aq_ex_name,sep=' ')
    df.to_csv("sfr_out.dat")

def setup():
    org_df = load_sfrout_aq_ex(sfr_file)
    org_df.to_csv(org_ex_name,sep=' ')
    segrech = list(org_df.index)
    kperkstp = list(org_df.columns)
    with open(aq_ex_name+".ins",'w') as f:
        f.write("pif ~\n")
        for i,sr in enumerate(segrech):
            if i == 0:
                f.write("l2 ")
            else:
                f.write("l1 ")
            for kk in kperkstp:
                obsnme = " w !" + sr + "_" + kk + '!'
                f.write(obsnme)
            f.write('\n')

if __name__ == "__main__":
    #setup()
    apply()