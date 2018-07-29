import os
import shutil
import numpy as np
import pandas as pd

import pyemu

MODEL_DIR = "Project_Default_ModelA"

# path to the model input directory
MODEL_INPUT_DIR = os.path.join(MODEL_DIR,"Inputs")

# path to the model output directory
MODEL_OUTPUT_DIR = os.path.join(MODEL_DIR,"Outputs")

# directories needed for multiplier process - not used
BASE_DIR = "base_inputs"
MULT_DIR = "mult_inputs"

# columns in reaches.csv to parameterize - no used, now using the par bound csv file
#PAR_COLS = ["RechargeFraction","GWDischargeFraction","SourceYield","GWMixDepth"]
PAR_BOUND_FILE = os.path.join("parameter_bound_set_up_reach.csv")
PB_DF = pd.read_csv(PAR_BOUND_FILE)
PB_DF.index = PB_DF.Parameter


# columns in reaches output csv to make observations with
#OBS_COLS = ["FlowOutTotal","FlowOutLateral","FluxOutTotalTperY","ConcentrationFW"]
OBS_COLS = ['ReachID', 'TotalRunoff', 'RechargeFraction', 'StreamLossFraction',
              'GWDischargeFraction', 'SourceYield','RechargeConcFactor',
              'GWDecayK', 'GWDeepDecayK', 'GWDecayFraction', 'GWDeepDecayFraction',
              'GWMixDepth','FlowOutLateral', 'FluxOutLateralTperY',
              'ConcentrationMedian', 'AgeY', 'AgeGWY', 'GWExportFraction']


def setup_org_dir():
    """copy the original model inputs for safe keeping"""
    if os.path.exists(MODEL_DIR):
        shutil.rmtree(MODEL_DIR)
    shutil.copytree(os.path.join(MODEL_DIR+"_org"),MODEL_DIR)




def setup_reach_tpl_file():
    """write a template file for
    the reaches.csv file - saves existing values to par.csv
    Mangles the names to fit in the 12-char limit"""

    base_df = pd.read_csv(os.path.join(MODEL_DIR+"_org","Inputs", "Reaches.csv"), index_col=0)
    pnames, pvals = [], []
    par_ub,par_lb = [],[]
    dist = []
    #for p in PAR_COLS:
    used = []
    for p in PB_DF.Parameter:
        assert p in base_df.columns, p
        bname = p
        if bname in used:
            i = 0
            while True:
                bname = bname[:-1] + str(i)
                if bname not in used:
                    break
                if i == 9:
                    raise Exception()
                i += 1
        used.append(bname)
        names = base_df.index.map(lambda x: "{0}{1}".format(bname,x))
        pvals.extend(list(base_df.loc[:,p].copy().values))
        if "rech" in p.lower():
            print(base_df.loc[:,p])
            print(pvals[-base_df.shape[0]:])
            print()
        pnames.extend(list(names.str.lower().values))

        mx = float(PB_DF.loc[p,"Max"])
        mn = float(PB_DF.loc[p,"Min"])
        if PB_DF.loc[p,"IsMultiplier"] == 1:
            ub = list((base_df.loc[:,p] * mx).values)
            lb = list((base_df.loc[:,p] * mn).values)
        else:
            ub = [PB_DF.loc[p,"Max"]] * names.shape[0]
            lb = [PB_DF.loc[p, "Min"]] * names.shape[0]
        par_ub.extend(ub)
        par_lb.extend(lb)
        dist.extend([PB_DF.loc[p,"Distribution"].lower().strip()] * names.shape[0])
        base_df.loc[:, p] = names.map(lambda x: "~  {0}  ~".format(x))

    with open("reaches.csv.tpl",'w') as f:
        f.write("ptf ~\n")
        f.flush()
        base_df.to_csv(f,mode='a')
    print(len(pnames))
    print(len(pvals))

    df = pd.DataFrame(data={"parnme":pnames,"parval1":pvals,
                            "parubnd":par_ub,"parlbnd":par_lb,
                            "dist":dist},index=pnames)

    df.to_csv("par.csv")


def setup_output_ins_files():
    """write the instruction file for outputs based on columns listed in OBS_COLS.
    mangles the names to fit in the 20-char limit"""
    f_out = open("processed_reach_output.dat",'r')
    f_ins = open("processed_reach_output.dat.ins",'w')
    f_ins.write("pif ~\n")
    header = f_out.readline().strip().split()
    used = []
    prefixes = {}
    for ih,h in enumerate(header[1:]):
        if h not in OBS_COLS:
            continue
        prefix = h.lower()
        if prefix in used:
            i = 1
            while True:
                prefix = prefix[:-1] + str(i)
                print(h[:10].lower(),prefix,used)
                if prefix not in used:
                    break
                if i == 9:
                    raise Exception(prefix)
                i += 1
        used.append(prefix)
        prefixes[ih] = prefix
    f_ins.write('l1\n')
    used = []
    for line in f_out:
        raw = line.strip().split()
        reach = raw[0]
        f_ins.write('l1 ')
        for ih,h in enumerate(header[1:]):
            #if "deep" in h or "point" in h or "additional" in h:
            #    f_ins.write(" w ")
            #    continue
            #prefix = 'sw'
            #if "gw" in h:
            #    prefix = 'gw'
            if h not in OBS_COLS:
                f_ins.write(" w ")
                continue
            prefix = prefixes[ih]
            oname = "{0}_{1}".format(prefix,reach)
            f_ins.write(" w !{0}! ".format(oname))
        f_ins.write('\n')
    f_out.close()
    f_ins.close()
    #os.system("inschek {0}.ins {0}".format("processed_reach_output.dat"))


def setup_control_file():
    """write the base pst file."""

    tpl_files = ["reaches.csv.tpl"]
    in_files = [os.path.join(MODEL_INPUT_DIR,"Reaches.csv")]
    ins_files = ["processed_reach_output.dat.ins"]
    out_files = ["processed_reach_output.dat"]

    par_df = pd.read_csv("par.csv", index_col=0)
    pst = pyemu.Pst.from_io_files(tpl_files,in_files,ins_files,out_files)
    pst.model_command = "python forward_run.py"
    par = pst.parameter_data
    par.loc[par_df.index,"parval1"] = par_df.parval1
    par.loc[par_df.index, "parubnd"] = par_df.parubnd
    par.loc[par_df.index, "parlbnd"] = par_df.parlbnd


    par.loc[:,"partrans"] = "none"
    par.loc[par.apply(lambda x: x.parubnd == x.parlbnd, axis=1), "partrans"] = "fixed"
    pst.control_data.noptmax = 0
    pst.write("clues.pst")
    pyemu.os_utils.run("pestpp clues.pst")


def draw():
    num_reals = 10000
    pst = pyemu.Pst(os.path.join("clues.pst"))
    par_df = pd.read_csv("par.csv", index_col=0)
    pe = pyemu.ParameterEnsemble.from_mixed_draws(pst=pst,
                                                  how_dict=par_df.dist.to_dict(),
                                                  num_reals=num_reals)
    #pe.to_csv("sweep_in.csv")
    pe.to_binary("sweep_in.jcb")

def invest():
    df = pd.read_csv("processed_reach_output.dat",delim_whitespace=True)
    print(df.shape)
    print(df.dropna().shape)
    print(df.head())
    print(df.dropna().head())

def invest2():
    df = pd.read_csv("par.csv",index_col=0)
    print(df.columns)
    print(df.loc[df.parnme.apply(lambda x: "rech" in x),"parval1"])

    df_org = pd.read_csv(os.path.join("..",MODEL_INPUT_DIR,"Reaches.csv"),index_col=0)
    df_new = pd.read_csv(os.path.join(MODEL_INPUT_DIR,"Reaches.csv"),index_col=0)
    print(df_org.shape,df_new.shape)
    diff = df_org - df_new
    print(diff.max())
    for col in diff.columns:
        if np.abs(diff.loc[:,col].max()) > 1.0e-6:
            idx = diff.loc[:,col].apply(np.abs) > 1.0e-6
            print(df_org.loc[idx,col])
            print(df_new.loc[idx,col])
            print()

if __name__ == "__main__":
    setup_org_dir()
    import forward_run
    setup_reach_tpl_file()
    setup_output_ins_files()
    setup_control_file()
    # invest2()
    draw()