# TODO: test variance and mean of draws, add chenoliver and test approx and full solution
import os
import shutil
import platform
import numpy as np
import pandas as pd
import pyemu

tests = """0) 10par_xsec "standard user mode" - draw reals from par-bounds prior and obs noise from weights
0a) 10par_xsec same as 0) but with multple lambda 
1) 10par_xsec start with existing par csv and obs csv - using empirical parcov and obscov
1a) 10par_xsec start with existing par csv and obs csv - using parcov file
2) 10par_xsec start with existing par csv and drawing obs en from weights 
3) 10par_xsec restart with full simulated obs en
3a) 10par_xsec restart with failed runs in simulated obs en
3b) 10par_xsec restart with failed runs and bad phi runs in simulated obs en with multiple lam
4) 10par_xsec reg_factor = 0.5 test
5)  10par_xsec full solution test with standard draw mode
5a) 10par_xsec full solution test with empirical parcov
6) freyberg "standard user mode" - draw reals from par-bounds prior and obs noise from weights
6a) freyberg same as 0) but with multple lambda 
7) freyberg draw par en from full parcov supplied in file
8) freyberg full solution with empirical parcov - supplied par csv, obs csv and restart csv with fails, bad phi,MAP solution, prior scaling, lam mults 
9) synth restart and upgrade 1.1M par problem"""

ies_vars = ["ies_par_csv", "ies_obs_csv", "ies_restart_obs_csv",
            "ies_bad_phi", "parcov_filename", "ies_num_reals",
            "ies_use_approx", "ies_use_prior_scaling", "ies_reg_factor",
            "ies_lambda_mults", "ies_init_lambda"]

exe_path = os.path.join("..", "..", "..", "exe", "windows", "x64", "Release", "pestpp-ies.exe")

noptmax = 3

compare_files = ["pest.phi.actual.csv", "pest.phi.meas.csv", "pest.phi.regul.csv",
                 "pest.{0}.par.csv".format(noptmax), "pest.{0}.obs.csv".format(noptmax)]
diff_tol = 1.0e-6


def write_empty_test_matrix():
    test_names = [t.split()[0] for t in tests.split('\n')]
    df = pd.DataFrame(index=test_names, columns=ies_vars)
    df.loc[:, "text"] = tests.split('\n')
    df.to_csv("ies_test.blank.csv")


def setup_suite_dir(model_d):
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))

    # set first par as fixed
    pst.parameter_data.loc[pst.par_names[0], "partrans"] = "fixed"

    # set noptmax
    pst.control_data.noptmax = noptmax

    # wipe all pestpp options
    pst.pestpp_options = {}

    # write a generic 2D cov
    if not os.path.exists(os.path.join(new_d, "prior.cov")):
        cov = pyemu.Cov.from_parameter_data(pst)
        cov = pyemu.Cov(cov.as_2d, names=cov.row_names)
        cov.to_ascii(os.path.join(new_d, "prior.cov"))
        cov.to_binary(os.path.join(new_d, "prior.jcb"))
    else:
        cov = pyemu.Cov.from_ascii(os.path.join(new_d, "prior.cov"))
    # draw some ensembles
    num_reals = 100
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=cov, num_reals=num_reals)
    pe.to_csv(os.path.join(new_d, "par.csv"))
    pe.to_binary(os.path.join(new_d, "par.jcb"))
    pe.to_csv(os.path.join(new_d, "sweep_in.csv"))
    pe.loc[:, pst.adj_par_names].to_csv(os.path.join(new_d, "par_some.csv"))

    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst, num_reals=num_reals)
    oe.to_csv(os.path.join(new_d, "obs.csv"))
    oe.to_binary(os.path.join(new_d, "obs.jcb"))

    pst.write(os.path.join(new_d, "pest.pst"))
    # run sweep
    if os.path.exists("master_sweep"):
        shutil.rmtree("master_sweep")
    pyemu.helpers.start_slaves(new_d, "sweep", "pest.pst", 10, master_dir="master_sweep")

    # process sweep output as restart csv and jcb
    df = pd.read_csv(os.path.join("master_sweep", "sweep_out.csv"))
    df.columns = [c.lower() for c in df.columns]
    df.to_csv(os.path.join(new_d, "restart.csv"))
    df = df.loc[:, pst.nnz_obs_names]
    df.loc[:, pst.nnz_obs_names].to_csv(os.path.join(new_d, "restart_some.csv"))
    df.iloc[:-3, :].to_csv(os.path.join(new_d, "restart_failed.csv"))

    # run pyemu ies
    pyemu_d = os.path.join(model_d, "master_pyemu")
    if os.path.exists(pyemu_d):
        shutil.rmtree(pyemu_d)
    shutil.copytree(new_d, pyemu_d)
    bdir = os.getcwd()
    os.chdir(pyemu_d)
    ies = pyemu.EnsembleSmoother("pest.pst", num_slaves=10, verbose="ies.log", slave_dir=os.path.join("..", "template"))
    ies.initialize(parensemble="par.csv", obsensemble="obs.csv")
    for i in range(pst.control_data.noptmax):
        ies.update()
    os.chdir(bdir)

    # pestpp_d = os.path.join(model_d,"master")
    # if os.path.exists(pestpp_d):
    #     shutil.rmtree(pestpp_d)
    # shutil.copytree(new_d,pestpp_d)


def run_suite(model_d):
    df = pd.read_csv("ies_test.csv")
    df = df.loc[df.text.apply(lambda x: model_d in x), :]
    df.fillna('', inplace=True)
    print(df)
    #model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")

    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))

    for i in range(df.shape[0]):
    #for i in range(3,4):
    
        test_vars = df.iloc[i, :].to_dict()
        test_name = test_vars["text"].split()[0].replace(")", '')
        print(test_vars["text"])
        pst.pestpp_options = {}
        for v in ies_vars:
            if pd.notnull(test_vars[v]):
                pst.pestpp_options[v] = test_vars[v]
        pst.write(os.path.join(template_d, "pest.pst"))
        test_d = os.path.join(model_d, "master_test_{0}".format(test_name))

        if os.path.exists(test_d):
            try:
                shutil.rmtree(test_d)
            except:
                print("error removing existing test_d: {0}".format(test_d))
                continue
        pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                                   master_dir=test_d, verbose=True, slave_root=model_d)


def compare_suite(model_d):
    base_d = os.path.join(model_d, "baseline_opt")
    test_ds = [d for d in os.listdir(model_d) if "master_test" in d]
    errors = []
    for test_d in test_ds:
        test_d = os.path.join(model_d, test_d)
        for compare_file in compare_files:
            if not os.path.exists(os.path.join(test_d, compare_file)):
                errors.append("missing compare file '{0}' in test_d '{1}'".format(compare_file, test_d))
            else:
                base_file = os.path.join(base_d, "{0}__{1}".
                                         format(os.path.split(test_d)[-1], compare_file))
                test_file = os.path.join(test_d, compare_file)
                try:
                    base_df = pd.read_csv(base_file,index_col=0)
                except Exception as e:
                    errors.append("error loading base_file {0}: {1}".format(base_file, str(e)))
                    continue
                try:
                    test_df = pd.read_csv(test_file,index_col=0)
                except Exception as e:
                    errors.append("error loading test_file {0}, {1}: {2}".format(test_d, base_file, str(e)))
                    continue
                try:
                    diff = (test_df - base_df).apply(np.abs)
                except Exception as e:
                    errors.append("error differencing base and test for '{0}':{1}".format(base_file, str(e)))
                max_diff = diff.max().max()
                if max_diff > diff_tol:
                    errors.append("max diff greater than diff tol for '{0}':{1}".format(base_file, max_diff))
    if len(errors) > 0:
        for e in errors:
            print("ERROR: ", e)
        raise Exception("errors in {0}: ".format(model_d) + '\n'.join(errors))


def test_10par_xsec():
    run_suite("ies_10par_xsec")
    compare_suite("ies_10par_xsec")


def test_freyberg():
    run_suite("ies_freyberg")
    compare_suite("ies_freyberg")


def rebase(model_d):
    """reset the "base" for the standard test suite"""
    # run_suite(model_d)
    base_d = os.path.join(model_d, "baseline_opt")
    if os.path.exists(base_d):
        shutil.rmtree(base_d)
    os.mkdir(base_d)

    # find test dirs
    print(os.listdir(model_d))
    test_ds = [d for d in os.listdir(model_d) if "master_test" in d]
    for test_d in test_ds:
        test_d = os.path.join(model_d, test_d)
        for compare_file in compare_files:
            if not os.path.exists(os.path.join(test_d, compare_file)):
                print("WARNING missing compare file:", test_d, compare_file)
            else:
                shutil.copy2(os.path.join(test_d, compare_file),
                             os.path.join(base_d, "{0}__{1}".
                                          format(os.path.split(test_d)[-1], compare_file)))

def tenpar_full_cov_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_full_cov_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = os.path.join(test_d,"pest.pst")
    pst = pyemu.Pst(pst_name)
    par = pst.parameter_data
    par.loc[:,"partrans"] = "fixed"
    par.loc[pst.par_names[:2],"partrans"] = "none"

    x = np.zeros((pst.npar_adj,pst.npar_adj)) + 0.25
    for i in range(pst.npar_adj):
        x[i,i] = 1.0
    cov = pyemu.Cov(x,names=pst.adj_par_names)
    cov.to_ascii(os.path.join(test_d,"prior.cov"))
    num_reals = 100000
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals=num_reals,use_homegrown=True)
    pe.enforce()

    pst.control_data.noptmax = 0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.write(pst_name)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)

    df = pd.read_csv(os.path.join(test_d,"pest.0.par.csv"),index_col=0)
    df.columns = [c.lower() for c in df.columns]

    p1,p2 = pst.adj_par_names
    pe_corr = pe.corr().loc[p1,p2]
    df_corr = df.corr().loc[p1,p2]
    diff = np.abs((pe_corr - df_corr)/pe_corr)
    assert diff < 0.05,"{0},{1},{2}".format(pe_corr,df_corr,diff)

    par.loc[pst.adj_par_names,"partrans"] = "log"
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals=num_reals, use_homegrown=True)
    pe.enforce()
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]

    p1, p2 = pst.adj_par_names
    pe_corr = pe.corr().loc[p1, p2]
    df_corr = df.corr().loc[p1, p2]
    diff = np.abs((pe_corr - df_corr) / pe_corr)
    assert diff < 0.05, "{0},{1},{2}".format(pe_corr, df_corr, diff)


def tenpar_subset_test():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_subset_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(base_d,test_d)
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.control_data.noptmax = 3

    # first without subset
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_lambda_mults"] = "1.0"
    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                               slave_root=model_d, master_dir=test_d)
    df_base = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_lambda_mults"] = "1.0"
    pst.pestpp_options["ies_subset_size"] = 15
    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                               slave_root=model_d, master_dir=test_d)
    df_sub = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)
    diff = (df_sub - df_base).apply(np.abs)
    print(diff.max())
    assert diff.max().max() == 0.0

def test_freyberg_full_cov():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    # if os.path.exists(test_d):
    #     shutil.rmtree(test_d)
    # shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.parameter_data.loc[:,"partrans"] = "log"
    pst.control_data.noptmax = 0
    pst.pestpp_options = {}
    num_reals = 500
    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.write(os.path.join(template_d,"pest.pst"))
    cov = pyemu.Cov.from_binary(os.path.join(template_d,"prior.jcb"))
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals,use_homegrown=True)
    pe.to_csv(os.path.join(template_d,"pyemu_pe.csv"))

    # pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
    #                            slave_root=model_d, master_dir=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    df.columns = [c.lower() for c in df.columns]
    diff_tol = 0.1
    pe = pe.apply(np.log10)
    for c in df.columns:
        if c not in pe.columns:
            continue

        m1, m2 = pe.loc[:,c].mean(), df.loc[:,c].mean()
        s1,s2 =  pe.loc[:,c].std(), df.loc[:,c].std()
        mdiff = np.abs((m1 - m2))
        sdiff = np.abs((s1 - s2))
        print(c,mdiff,sdiff)
        assert mdiff < diff_tol,"mean fail {0}:{1},{2},{3}".format(c,m1,m2,mdiff)
        assert sdiff < diff_tol,"std fail {0}:{1},{2},{3}".format(c,s1,s2,sdiff)

    #look for bias
    diff = df - pe
    assert diff.mean().mean() < 0.01

def invest():
    d = os.path.join("ies_freyberg","master_draw_test")
    df = pd.read_csv(os.path.join(d,"draws.dat"),delim_whitespace=True,header=None)
    print(df.std().mean(),df.mean().mean())

    df1 = pd.read_csv(os.path.join(d,"pest.0.par.csv"), index_col=0)
    df2 = pd.read_csv(os.path.join(d,"pyemu_pe.csv"), index_col=0)
    df1.columns = [c.lower() for c in df1.columns]

    df1 = df1.apply(np.log10)
    df2 = df2.apply(np.log10)

    for p in df1.columns:
        print(p)
        #print(p)
        #print(df1.loc[:,p])
        #print(df2.loc[:,p])
        print(p,df1.loc[:, p].std(), df2.loc[:, p].std())
        #break

if __name__ == "__main__":
    # write_empty_test_matrix()
    #setup_suite_dir("ies_freyberg")
    #run_suite("ies_freyberg")
    #run_suite("ies_10par_xsec")
    #rebase("ies_freyberg")
    #rebase("ies_10par_xsec")
    #test_10par_xsec()
    #test_freyberg()
    #test_freyberg_full_cov()
    #invest()
    #compare_suite("ies_10par_xsec")
    #compare_suite("ies_freyberg")
    #tenpar_subset_test()
    tenpar_full_cov_test()
