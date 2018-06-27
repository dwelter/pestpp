# TODO: test variance and mean of draws, add chenoliver and test approx and full solution
import os
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import matplotlib.pyplot as plt
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
            "ies_lambda_mults", "ies_init_lambda","ies_include_base","ies_subset_size"]

if "windows" in platform.platform().lower():
    exe_path = os.path.join("..", "..", "..", "exe", "windows", "x64", "Release", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower():
    exe_path = os.path.join("..", "..", "..", "exe", "mac", "pestpp-ies")
else:
    exe_path = os.path.join("..", "..", "..", "exe", "linux", "pestpp-ies")

noptmax = 3

compare_files = ["pest.phi.actual.csv", "pest.phi.meas.csv", "pest.phi.regul.csv",
                 "pest.{0}.par.csv".format(noptmax), "pest.{0}.obs.csv".format(noptmax),
                 "pest.{0}.par.csv".format(0), "pest.base.obs.csv"]
diff_tol = 1.0e-6

num_reals = 30
def write_empty_test_matrix():
    test_names = [t.split()[0] for t in tests.split('\n')]
    df = pd.DataFrame(index=test_names, columns=ies_vars)
    df.loc[:, "text"] = tests.split('\n')
    df.to_csv("ies_test.blank.csv")


def setup_suite_dir(model_d):
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))

    # set first par as fixed
    #pst.parameter_data.loc[pst.par_names[0], "partrans"] = "fixed"

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = 1.0

    # set noptmax
    pst.control_data.noptmax = noptmax

    # wipe all pestpp options
    pst.pestpp_options = {}

    # write a generic 2D cov
    if os.path.exists(os.path.join(new_d,"prior.jcb")):
        cov = pyemu.Cov.from_binary(os.path.join(new_d,"prior.jcb"))
        cov.to_ascii(os.path.join(new_d,"prior.cov"))
    elif os.path.exists(os.path.join(new_d, "prior.cov")):
        cov = pyemu.Cov.from_ascii(os.path.join(new_d, "prior.cov"))
    else:
        cov = pyemu.Cov.from_parameter_data(pst)
        cov = pyemu.Cov(cov.as_2d, names=cov.row_names)
        cov.to_ascii(os.path.join(new_d, "prior.cov"))
        cov.to_binary(os.path.join(new_d, "prior.jcb"))

    # draw some ensembles
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=cov, num_reals=num_reals,
                                                    use_homegrown=True,group_chunks=True)
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
    pyemu.helpers.start_slaves(new_d, "pestpp-swp", "pest.pst", 10, master_dir="master_sweep",port=4020)

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
        # if test_vars["pyemu_compare"] == 0:
        #     continue
        test_name = test_vars["text"].split()[0].replace(")", '')
        print(test_vars["text"])
        #if "6" not in test_vars["text"]:
        #    continue
        pst.pestpp_options = {}
        for v in ies_vars:
            if pd.notnull(test_vars[v]):
                try:
                    pst.pestpp_options[v] = test_vars[v].replace('"','')
                except:
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
                                   master_dir=test_d, verbose=True, slave_root=model_d,
                                   port=4020)


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


def compare_pyemu():
    df = pd.read_csv("ies_test.csv")
   # df = df.loc[df.text.apply(lambda x: model_d in x), :]
    df = df.loc[df.pyemu_compare==1,:]
    df.fillna('', inplace=True)
    for i in range(df.shape[0]):
        args = df.iloc[i,:].to_dict()
        model_d = args["text"].split()[1]
        test_d = "master_test_"+args["text"].split()[0].replace(')','')
        print(model_d)
        pyemu_d = os.path.join(model_d,"master_pyemu")
        assert os.path.exists(pyemu_d),pyemu_d
        test_d = os.path.join(model_d,test_d)
        assert os.path.exists(test_d),test_d
        print(test_d)
        pp_phi_file = os.path.join(test_d,"pest.phi.actual.csv")
        py_phi_file = os.path.join(pyemu_d,"pest.pst.iobj.actual.csv")
        try:
            pp_phi = pd.read_csv(pp_phi_file)
        except Exception as e:
            raise Exception("error loading pp_phi {0}: {1}".format(pp_phi_file,str(e)))
        try:
            py_phi = pd.read_csv(py_phi_file)
        except Exception as e:
            raise Exception("error loading py_phi {0}: {1}".format(py_phi_file, str(e)))

        m_diff = pp_phi.loc[:,"mean"] - py_phi.loc[:,"mean"]
        print(m_diff)
        m_diff.dropna(inplace=True)
        assert m_diff.shape[0] == noptmax + 1
        assert m_diff.apply(np.abs).max() < 0.035
        s_diff = pp_phi.loc[:,"standard_deviation"] - py_phi.loc[:,"std"]
        print(s_diff)
        s_diff.dropna(inplace=True)
        assert s_diff.shape[0] == noptmax + 1
        assert s_diff.apply(np.abs).max() < 0.035,s_diff



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


def tenpar_narrow_range_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_narrow_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst_name = os.path.join(test_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    par = pst.parameter_data
    #par.loc[:, "partrans"] = "fixed"
    par.loc[:, "parubnd"] = 1.0e+10 #par.parval1 * 1.0001
    par.loc[:, "parlbnd"] = 1.0e-10 #par.parval1 * 0.9999
    #par.loc[pst.par_names[:2], "partrans"] = "none"
    #par.loc[pst.par_names[0],"pargp"] = "stage"

    x = np.zeros((pst.npar_adj, pst.npar_adj)) + 1.0e-11
    for i in range(pst.npar_adj):
        x[i, i] = 5.0e-10
    cov = pyemu.Cov(x, names=pst.adj_par_names)
    cov.to_ascii(os.path.join(test_d, "prior.cov"))
    num_reals = 5000
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals=num_reals, use_homegrown=True)
    #pe.enforce()
    #pe.to_csv(os.path.join(test_d,"pyemu_draws.csv"))

    pst.control_data.noptmax = 0

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] *= 1.5

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_group_draws"] = False
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)

    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]
    p1, p2 = pst.adj_par_names[:2]
    v1,v2 = pe.loc[:,p1].var(),df.loc[:,p1].var()
    diff = np.abs(100 * ((v1 - v2) / v1))
    print(v1,v2,diff)
    assert diff < 1.0

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_group_draws"] = True
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)

    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]
    p1, p2 = pst.adj_par_names[:2]
    v1, v2 = pe.loc[:, p1].var(), df.loc[:, p1].var()
    diff = np.abs(100 * ((v1 - v2) / v1))
    print(v1, v2, diff)
    assert diff < 1.0

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
    par.loc[:,"parubnd"] = 1.0e+10
    par.loc[:,"parlbnd"] = -1.0e+10
    par.loc[pst.par_names[:2],"partrans"] = "none"

    x = np.zeros((pst.npar_adj,pst.npar_adj)) + 0.25
    for i in range(pst.npar_adj):
        x[i,i] = 1.0
    cov = pyemu.Cov(x,names=pst.adj_par_names)
    cov.to_ascii(os.path.join(test_d,"prior.cov"))
    num_reals = 10000
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
    par.loc[:,"parlbnd"] = 1.0e-10
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals=num_reals, use_homegrown=True)
    pe.enforce()
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]

    p1, p2 = pst.adj_par_names
    pe_corr = pe.apply(np.log10).corr().loc[p1, p2]
    df_corr = df.apply(np.log10).corr().loc[p1, p2]
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
    pst.pestpp_options["ies_num_reals"] = 20
    pst.pestpp_options["ies_lambda_mults"] = "1.0"
    pst.pestpp_options["ies_accept_phi_fac"] = 100.0
    pst.pestpp_options["ies_subset_size"] = 21
    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                               slave_root=model_d, master_dir=test_d,port=4020)
    df_base = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 20
    pst.pestpp_options["ies_lambda_mults"] = "1.0"
    pst.pestpp_options["ies_subset_size"] = 5
    pst.pestpp_options["ies_accept_phi_fac"] = 100.0

    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                               slave_root=model_d, master_dir=test_d,port=4020)
    df_sub = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)
    diff = (df_sub - df_base).apply(np.abs)
    diff = diff.iloc[:,6:]
    print(diff.max())
    print(df_sub.iloc[-1,:])
    print(df_base.iloc[-1,:])
    assert diff.max().max() == 0.0

def test_freyberg_full_cov():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    pst.parameter_data.loc[:,"partrans"] = "log"
    
    pst.control_data.noptmax = 0
    pst.pestpp_options = {}
    num_reals = 5000

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"

    pst.write(os.path.join(test_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals, use_homegrown=True)
    pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    # pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
    #                            slave_root=model_d, master_dir=test_d)
    # pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    # print("loading df")
    # df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    # df.columns = [c.lower() for c in df.columns]
    # pe = pe.apply(np.log10)
    # pe_corr = pe.corr()
    # df_corr = df.corr()

    # diff_tol = 0.05

    # for c in df.columns:
    #     if c not in pe.columns:
    #         continue

    #     m1, m2 = pe.loc[:, c].mean(), df.loc[:, c].mean()
    #     s1, s2 = pe.loc[:, c].std(), df.loc[:, c].std()
    #     mdiff = np.abs((m1 - m2))
    #     sdiff = np.abs((s1 - s2))
    #     print(c, mdiff, sdiff)
    #     assert mdiff < diff_tol, "mean fail {0}:{1},{2},{3}".format(c, m1, m2, mdiff)
    #     assert sdiff < diff_tol, "std fail {0}:{1},{2},{3}".format(c, s1, s2, sdiff)

    # # look for bias
    # diff = df - pe
    # assert diff.mean().mean() < 0.01


    #full cov
    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.pestpp_options["ies_group_draws"] = 'true'
    pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"
    pst.write(os.path.join(test_d,"pest.pst"))
    cov = pyemu.Cov.from_binary(os.path.join(test_d,"prior.jcb"))

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals,use_homegrown=True)
    pe.to_csv(os.path.join(test_d,"pyemu_pe.csv"))

    # pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
    #                            slave_root=model_d, master_dir=test_d)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    df.columns = [c.lower() for c in df.columns]
    pe = pe.apply(np.log10)
    pe_corr = pe.corr()
    df_corr = df.corr()

    for i,p1 in enumerate(pst.adj_par_names):
        for p2 in pst.adj_par_names[i+1:]:
            c1 = pe_corr.loc[p1,p2]
            c2 = df_corr.loc[p1,p2]
            #print(p1,p2,c1,c2)

    diff_tol = 0.05

    for c in df.columns:
        if c not in pe.columns:
            continue

        m1, m2 = pe.loc[:,c].mean(), df.loc[:,c].mean()
        s1,s2 =  pe.loc[:,c].std(), df.loc[:,c].std()
        mdiff = np.abs((m1 - m2))
        sdiff = np.abs((s1 - s2))
        #print(c,mdiff,sdiff)
        assert mdiff < diff_tol,"mean fail {0}:{1},{2},{3}".format(c,m1,m2,mdiff)
        assert sdiff < diff_tol,"std fail {0}:{1},{2},{3}".format(c,s1,s2,sdiff)

    #look for bias
    diff = df - pe
    assert diff.mean().mean() < 0.01


def test_freyberg_full_cov_reorder():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    par.loc[pst.par_names[:5],"pargp"] = pst.par_groups[-1]

    pst.control_data.noptmax = 0
    pst.pestpp_options = {}
    num_reals = 1000

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"

    pst.write(os.path.join(test_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals, use_homegrown=True)
    pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    # pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
    #                            slave_root=model_d, master_dir=test_d)
    # pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    # print("loading df")
    # df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    # df.columns = [c.lower() for c in df.columns]
    # pe = pe.apply(np.log10)
    # pe_corr = pe.corr()
    # df_corr = df.corr()

    # diff_tol = 0.05

    # for c in df.columns:
    #     if c not in pe.columns:
    #         continue

    #     m1, m2 = pe.loc[:, c].mean(), df.loc[:, c].mean()
    #     s1, s2 = pe.loc[:, c].std(), df.loc[:, c].std()
    #     mdiff = np.abs((m1 - m2))
    #     sdiff = np.abs((s1 - s2))
    #     print(c, mdiff, sdiff)
    #     assert mdiff < diff_tol, "mean fail {0}:{1},{2},{3}".format(c, m1, m2, mdiff)
    #     assert sdiff < diff_tol, "std fail {0}:{1},{2},{3}".format(c, s1, s2, sdiff)

    # # look for bias
    # diff = df - pe
    # assert diff.mean().mean() < 0.01


    #full cov
    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.pestpp_options["ies_group_draws"] = 'true'
    pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"
    pst.write(os.path.join(test_d,"pest.pst"))
    cov = pyemu.Cov.from_binary(os.path.join(test_d,"prior.jcb"))

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals,use_homegrown=True)
    pe.to_csv(os.path.join(test_d,"pyemu_pe.csv"))

    # pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
    #                            slave_root=model_d, master_dir=test_d)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    df.columns = [c.lower() for c in df.columns]
    pe = pe.apply(np.log10)
    pe_corr = pe.corr()
    df_corr = df.corr()

    for i,p1 in enumerate(pst.adj_par_names):
        for p2 in pst.adj_par_names[i+1:]:
            c1 = pe_corr.loc[p1,p2]
            c2 = df_corr.loc[p1,p2]
            #print(p1,p2,c1,c2)

    diff_tol = 0.05

    for c in df.columns:
        if c not in pe.columns:
            continue

        m1, m2 = pe.loc[:,c].mean(), df.loc[:,c].mean()
        s1,s2 =  pe.loc[:,c].std(), df.loc[:,c].std()
        mdiff = np.abs((m1 - m2))
        sdiff = np.abs((s1 - s2))
        #print(c,mdiff,sdiff)
        assert mdiff < diff_tol,"mean fail {0}:{1},{2},{3}".format(c,m1,m2,mdiff)
        assert sdiff < diff_tol,"std fail {0}:{1},{2},{3}".format(c,s1,s2,sdiff)

    #look for bias
    diff = df - pe
    assert diff.mean().mean() < 0.01


def test_freyberg_full_cov_reorder_run():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    #par.loc[pst.par_names[:5],"pargp"] = pst.par_groups[-1]

    pst.control_data.noptmax = 10
    pst.pestpp_options = {}
    num_reals = 100

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"

    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "true"
    pst.pestpp_options["ies_group_draws"] = 'true'
   # pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"
    
    pst.write(os.path.join(template_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    #pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals, use_homegrown=True)
    #pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=10,
                                slave_root=model_d, master_dir=test_d,port=4020)
    

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

def test_synth():
    model_d = "ies_synth"
    test_d = os.path.join(model_d,"master")
    template_d = os.path.join(model_d,"template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.pestpp_options = {}
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.9,1.1]
    pst.pestpp_options["ies_num_reals"] = 30
    pst.pestpp_options["ies_save_binary"] = True
    pst.control_data.noptmax = 2
    print("writing pst")
    pst.write(os.path.join(template_d,"pest.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d,exe_path,"pest.pst",num_slaves=15,
        master_dir=test_d,slave_root=model_d,port=4020)

def test_chenoliver():
    model_d = "ies_chenoliver"
    test_d = os.path.join(model_d,"master")
    template_d = os.path.join(model_d,"template")

    # build the prior cov matrix
    x = np.zeros((1,1)) + 0.5
    cov = pyemu.Cov(x=x,names=["par"])
    cov.to_ascii(os.path.join(template_d,"prior.cov"))

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.parameter_data.loc[:,"parlbnd"] = -1.0e+10
    pst.parameter_data.loc[:,"parubnd"] = 1.0+10
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = 100
    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    
    num_reals = 500
    noptmax = 4
    

    shutil.rmtree(test_d)
    
    pst.observation_data.loc[:,"weight"] = 0.25

    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax

    pst.write(os.path.join(template_d,"pest.pst"))
    

    pyemu.helpers.start_slaves(template_d,exe_path,"pest.pst",num_slaves=20,
        master_dir=test_d,slave_root=model_d,port=4020,silent_master=True)
    df_full_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_full_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.helpers.start_slaves(template_d,exe_path,"pest.pst",num_slaves=20,
        master_dir=test_d,slave_root=model_d,port=4020,silent_master=True)
    df_approx_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_approx_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    ax = plt.subplot(211)
    ax2 = plt.subplot(212)
    #ax.plot(df_full.loc[:,"mean"],color='b',label="full")
    #ax.plot(df_approx.loc[:,"mean"],color='g',label="approx")
    df_full_obs.OBS.hist(bins=30,color='b',alpha=0.5, ax=ax)
    df_approx_obs.OBS.hist(bins=30,color='0.5',alpha=0.5,ax=ax)
    df_full_par.PAR.hist(bins=30,color='b',alpha=0.5, ax=ax2)
    df_approx_par.PAR.hist(bins=30,color='0.5',alpha=0.5,ax=ax2)
    ax.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_obs.OBS.mean(),df_approx_obs.OBS.mean()))
    ax2.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_par.PAR.mean(),df_approx_par.PAR.mean()))
    plt.tight_layout()
    plt.savefig(os.path.join(model_d,"full_approx_ov16.png"))
    plt.close("all")
    # d = np.abs(df_full_par.PAR.mean() - 5.8)
    # assert d < 0.05,d
    # d = np.abs(df_full_par.PAR.mean() - 6.0)
    # assert d < 0.05,d

    pst.observation_data.loc[:,"weight"] = 1.0

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.helpers.start_slaves(template_d,exe_path,"pest.pst",num_slaves=25,
        master_dir=test_d,slave_root=model_d,port=4020,silent_master=True)
    df_full_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_full_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.helpers.start_slaves(template_d,exe_path,"pest.pst",num_slaves=25,
        master_dir=test_d,slave_root=model_d,port=4020,silent_master=True)
    df_approx_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_approx_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    ax = plt.subplot(211)
    ax2 = plt.subplot(212)
    #ax.plot(df_full.loc[:,"mean"],color='b',label="full")
    #ax.plot(df_approx.loc[:,"mean"],color='g',label="approx")
    df_full_obs.OBS.hist(bins=30,color='b',alpha=0.5, ax=ax)
    df_approx_obs.OBS.hist(bins=30,color='0.5',alpha=0.5,ax=ax)
    df_full_par.PAR.hist(bins=30,color='b',alpha=0.5, ax=ax2)
    df_approx_par.PAR.hist(bins=30,color='0.5',alpha=0.5,ax=ax2)
    ax.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_obs.OBS.mean(),df_approx_obs.OBS.mean()))
    ax2.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_par.PAR.mean(),df_approx_par.PAR.mean()))
    plt.tight_layout()
    plt.savefig(os.path.join(model_d,"full_approx_ov1.png"))
    plt.close("all")

    # d = np.abs(df_full_par.PAR.mean() - 5.99)
    # assert d < 0.05,d
    # d = np.abs(df_full_par.PAR.mean() - 6.0)
    # assert d < 0.05,d






def test_kirishima():

    model_d = "ies_kirishima"
    test_d = os.path.join(model_d, "master")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.pestpp_options = {}
    #pst.pestpp_options["ies_num_reals"] = 300
    #pst.pestpp_options["ies_use_approx"] = "true"
    #pst.pestpp_options["ies_use_prior_scaling"] = "true"
    #pst.pestpp_options["ies_subset_size"] = 10
    #pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]

    #pst.pestpp_options["ies_initial_lambda"] = 1000.0
    #pst.pestpp_options["ies_bad_phi"] = 80000.0
    pst.control_data.noptmax = 10
    print("writing pst")
    pst.write(os.path.join(template_d, "pest.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=15, master_dir=test_d,
                               slave_root=model_d,port=4020)

def test_freyberg_ineq():

    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "test_ineq")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.observation_data.loc[pst.nnz_obs_names[3],"obgnme"] = "less_than"
    pst.observation_data.loc[pst.nnz_obs_names[3],"weight"] = 100.0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 100
    pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]
    pst.control_data.noptmax = 3
    print("writing pst")
    pst.write(os.path.join(template_d, "pest_ineq.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_ineq.pst", num_slaves=10, master_dir=test_d,
                               slave_root=model_d,port=4020)

    obs_csvs = [f for f in os.listdir(test_d) if "obs" in f and f.endswith(".csv")]
    print(obs_csvs)
    df = pd.read_csv(os.path.join(test_d,obs_csvs[-1]),index_col=0)
    df.columns = df.columns.map(str.lower)
    pst = pyemu.Pst(os.path.join(template_d, "pest_ineq.pst"))

    axes = [plt.subplot(pst.nnz_obs,1,i+1) for i in range(pst.nnz_obs)]
    #df.loc[:,pst.nnz_obs_names].hist(bins=10,axes=axes)
    for i,n in enumerate(pst.nnz_obs_names):
        v = pst.observation_data.loc[n,"obsval"]
        ax = axes[i]
        df.loc[:,n].hist(ax=ax,bins=10)
        ax.plot([v,v],ax.get_ylim(),"k--",lw=2.0)
    #plt.show()

def tenpar_fixed_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "test_fixed")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.control_data.noptmax = 2
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=50)

    pe.loc[:,"stage"] = np.linspace(0.0,1.0,pe.shape[0])
    pe.loc[:,"k_01"] = 5.0
    pe.to_csv(os.path.join(template_d,"par.csv"))
    pe.to_binary(os.path.join(template_d, "par.jcb"))
    fixed_pars = ["stage","k_01"]
    pst.parameter_data.loc[fixed_pars,"partrans"] = "fixed"

    def compare():
        csvs = [f for f in os.listdir(test_d) if f.endswith(".par.csv")]
        dfs = [pd.read_csv(os.path.join(test_d,csv),index_col=0) for csv in csvs]
        for df in dfs:
            df.columns = df.columns.map(str.lower)
            df = df.loc[:,fixed_pars]
            df = df.iloc[:-1,:]
            df.index = pe.index
            diff = pe.loc[df.index,fixed_pars] - df
            assert diff.apply(np.abs).sum().sum() < 0.01, diff

    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.write(os.path.join(template_d, "pest.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=20, master_dir=test_d,
                               slave_root=model_d,port=4020)
    compare()
    pe.to_binary(os.path.join(template_d,"par.jcb"))
    pst.pestpp_options["ies_par_en"] = "par.jcb"
    pst.write(os.path.join(template_d, "pest.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=20, master_dir=test_d,
                               slave_root=model_d,port=4020)
    compare()

    pst.pestpp_options["ies_par_en"] = "par.jcb"
    pst.pestpp_options["ies_save_binary"] = 'true'
    pst.write(os.path.join(template_d, "pest.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest.pst", num_slaves=20, master_dir=test_d,
                               slave_root=model_d,port=4020)
    pe1 = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest.0.par.jcb")).iloc[:-1,:]
    pe1.index = pe.index
    diff = pe - pe1
    assert diff.apply(np.abs).sum().sum() == 0.0


def tenpar_weights_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "test_weights")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)

    dfs = []
    

    for i in range(3):
        obs = pst.observation_data.weight.copy()    
        dfs.append(obs)

    df = pd.concat(dfs,axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(test_d,"weights.csv"))
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst=pst,num_reals=df.shape[0])
    oe.to_csv(os.path.join(test_d,"obs.csv"))
    pst.control_data.noptmax = -1
    #pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.pestpp_options["ies_num_reals"] = df.shape[0]
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_obs_en"] = "obs.csv"

    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    df_act = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"))
    df_meas = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))

    pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    df_act1 = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"))
    df_meas1 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))
    print(df_act.loc[0,"mean"],df_act1.loc[0,"mean"])
    assert df_act.loc[0,"mean"] == df_act1.loc[0,"mean"]

    print(df_meas.loc[0,"mean"],df_meas1.loc[0,"mean"])
    assert df_meas.loc[0,"mean"] == df_meas1.loc[0,"mean"]


def tenpar_tight_tol_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "test_tighttol")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)

    
    pst.control_data.noptmax = 3
    #pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 0.00000001
    pst.pestpp_options["lambda_scale_fac"] = 10.0
    pst.pestpp_options["ies_initial_lambda"] = 0.000001
    pst.pestpp_options["ies_accept_phi_fac"] = 1.0

    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    

def tenpar_weight_pareto_test():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "test_weight_pareto")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
    #if os.path.exists(test_d):
    #   shutil.rmtree(test_d)
    #shutil.copytree(template_d,test_d)

    dfs = []
    obs = pst.observation_data
    obs.loc["h02_08", "obsval"] = 0
    obs.loc["h02_08", "weight"] = 1.0
    dfs = []
    weights = np.linspace(0.0,10.0,50)
    for weight in weights:
        obs = pst.observation_data.weight.copy()
        obs.loc["h02_08"] = weight
        dfs.append(obs)
        # for i in range(nreal_per):
        #     dfs.append(obs)
        #for weight2 in weights:
        #    obs1 = obs.copy()
        #    obs1.loc["h02_08"] = weight2
        #    #print(obs1)
        #    dfs.append(obs1)

    df = pd.concat(dfs,axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(template_d,"weights.csv"))

    obs = pst.observation_data.obsval
    df = pd.concat([obs.copy() for i in range(df.shape[0])],axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(template_d, "obs.csv"))


    pst.control_data.noptmax = 3
    pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.pestpp_options["ies_obs_en"] = "obs.csv"
    pst.pestpp_options["ies_num_reals"] = df.shape[0]
    pst.pestpp_options["ies_subset_size"] = df.shape[0]
    pst.pestpp_options["ies_save_binary"]= False
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    obs = pst.observation_data
    #obs.loc[pst.nnz_obs_names,"obsval"] = obs.loc[pst.nnz_obs_names,"obsval"].mean()
    #obs.loc["h01_04","obsval"] = 0.0
    #obs.loc["h01_06","obsval"] = 6.0

    pst.write(os.path.join(template_d,"pest_pareto.pst"))
    pyemu.os_utils.start_slaves(template_d,exe_path,"pest_pareto.pst",num_slaves=40,
                                slave_root=model_d,master_dir=test_d,port=4020)
    obs = pst.observation_data
    df_init = pd.read_csv(os.path.join(test_d,"pest_pareto.0.obs.csv".format(pst.control_data.noptmax)))
    df = pd.read_csv(os.path.join(test_d,"pest_pareto.{0}.obs.csv".format(pst.control_data.noptmax)))
    df_phi = pd.read_csv(os.path.join(test_d,"pest_pareto.phi.meas.csv"))
    df_phi = df_phi.loc[:,[str(v) for v in df_init.index.values]]
    # fig = plt.figure(figsize=(10,10))
    # ax = plt.subplot(221)
    # ax2 = plt.subplot(223)
    # ax3 = plt.subplot(224)
    # ax2.scatter(((df_init.H01_04-obs.loc["h01_04","obsval"])**2),
    #     ((df_init.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='0.5')
    # ax2.scatter(((df.H01_04-obs.loc["h01_04","obsval"])**2),
    #     ((df.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='b')
    # ax.scatter(df_phi.iloc[0,:],((df_init.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1,:],((df.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='b')
    # ax3.scatter(df_phi.iloc[0,:],((df_init.H01_04-obs.loc["h01_04","obsval"])**2),s=4, color='0.5')
    # ax3.scatter(df_phi.iloc[-1,:],((df.H01_04-obs.loc["h01_04","obsval"])**2),s=4, color='b')

    # ax2.scatter(df_init.H01_04,df_init.H01_06, s=4, color='0.5')
    # ax2.scatter(df.H01_04,df.H01_06, s=4, color='b')
    # ax.scatter(df_phi.iloc[0, :], df_init.H01_06, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], df.H01_06, s=4, color='b')
    # ax3.scatter(df_phi.iloc[0, :], df_init.H01_04, s=4, color='0.5')
    # ax3.scatter(df_phi.iloc[-1, :], df.H01_04, s=4, color='b')

    # fig = plt.figure(figsize=(10, 10))
    # ax = plt.subplot(221)
    # ax.scatter(df_phi.iloc[0, :], df_init.H02_08, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], df.H02_08, s=4, color='b')
    # ax.scatter(df_phi.iloc[0, :], (df_init.H02_08-obs.obsval["h02_08"])**2, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], (df.H02_08-obs.obsval["h02_08"])**2, s=4, color='b')

    #plt.show()


def rosenbrock_function():
    par_df = pd.read_csv("par.dat",delim_whitespace=True,index_col=0)
    tot = 0.0
    for i in range(par_df.shape[0]-1):
        tot += 100.0*(par_df.iloc[i+1] - par_df.iloc[i]**2)**2 + (1 - par_df.iloc[i])**2
    with open("obs.dat",'w') as f:
        f.write("obs {0:20.8E}".format(tot))

def setup_rosenbrock():

    npar = 2
    test_d = "ies_rosenbrock"
    if not os.path.exists(test_d):
        os.mkdir(test_d)

    template_d = os.path.join(test_d,"template")
    if os.path.exists(template_d):
        shutil.rmtree(template_d)
    os.mkdir(template_d)

    with open(os.path.join(template_d,"par.dat.tpl"),'w') as f:
        f.write("ptf ~\n")
        for i in range(npar):
            f.write("par{0:04d}  ~   par{0:04d}   ~\n")
    with open(os.path.join(template_d,"obs.dat.ins"),'w') as f:
        f.write("pif ~\n")
        f.write("l1 w !obs1!\n")

    bd = os.getcwd()




def tenpar_localizer_test1():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test1")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    #shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    
    #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat = pyemu.Matrix.from_names(["head"],pst.adj_par_names).to_dataframe()
    mat.loc[:,:] = 1.0
    #mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.control_data.noptmax = 3
    
    #pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d,"pest_local.pst")
    pst.write(pst_name)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11,
                                   master_dir=test_d, verbose=True, slave_root=model_d,
                                   port=4020)
    phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local.phi.actual.csv"))

    pst.pestpp_options.pop("ies_localizer")
    pst.write(pst_name)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11,
                                   master_dir=test_d, verbose=True, slave_root=model_d,
                                   port=4020)
    phi_df2 = pd.read_csv(os.path.join(test_d,"pest_local.phi.actual.csv"))
    diff = phi_df1 - phi_df2
    print(diff.max().max())
    assert diff.max().max() == 0
    plt.plot(phi_df1.total_runs,phi_df1.loc[:,"mean"], label="local")
    plt.plot(phi_df2.total_runs,phi_df2.loc[:,"mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d,"local_test.pdf"))

def tenpar_localizer_test2():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test2")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    #shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d,"par.csv"))

    #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    use_pars = pst.adj_par_names[0:2]

    #mat = pyemu.Matrix.from_names(["head"],use_pars).to_dataframe()
    mat = pyemu.Matrix.from_names(["head"],pst.adj_par_names).to_dataframe()
    
    mat.loc[:,:] = 0.0
    mat.loc[:,use_pars] = 1.0
    #mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_include_base"] = False
    pst.control_data.noptmax = 3
    
    #pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d,"pest_local.pst")
    pst.write(pst_name)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11,
                                   master_dir=test_d, verbose=True, slave_root=model_d,
                                   port=4020)
    phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local.phi.actual.csv"))
    par_df1 = pd.read_csv(os.path.join(test_d,"pest_local.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    par_df1.columns = par_df1.columns.str.lower()

    pst.pestpp_options.pop("ies_localizer")
    pst.parameter_data.loc[:,"partrans"] = "fixed"
    pst.parameter_data.loc[use_pars,"partrans"] = "log"
    pst.write(os.path.join(template_d,"pest_base.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_base.pst", num_slaves=11,
                                   master_dir=test_d+"_base", verbose=True, slave_root=model_d,
                                   port=4020)
    phi_df2 = pd.read_csv(os.path.join(test_d+"_base","pest_base.phi.actual.csv"))
    par_df2 = pd.read_csv(os.path.join(test_d+"_base","pest_base.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    par_df2.columns = par_df1.columns.str.lower()
    plt.plot(phi_df1.total_runs,phi_df1.loc[:,"mean"], label="local")
    plt.plot(phi_df2.total_runs,phi_df2.loc[:,"mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d,"local_test.pdf"))

    for par in par_df1.columns:
        diff = par_df1.loc[:,par] - par_df2.loc[:,par]
        print(diff)
        print(diff.sum())
        diff = par_df1.loc[:,par] - pe.loc[:,par]
        print(diff.sum())
        assert diff.sum() == 0.0
    diff = phi_df1 - phi_df2
    print(diff.max().max())
    assert diff.max().max() == 0


def tenpar_localizer_test3():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test3")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par.csv"))

    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    use_pars = pst.adj_par_names[0:2]


    mat.loc[:, :] = 0.0
    mat.loc[pst.nnz_obs_names, use_pars] = 1.0
    # mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_accept_phi_fac"] = 1.01
    pst.control_data.nphistp = 10
    pst.control_data.nphinored = 10
    pst.control_data.noptmax = 10

    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_local.pst")
    pst.write(pst_name)
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11,
                               master_dir=test_d, verbose=True, slave_root=model_d,
                               port=4020)
    phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))
    par_df1 = pd.read_csv(os.path.join(test_d, "pest_local.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df1.columns = par_df1.columns.str.lower()

    pst.pestpp_options.pop("ies_localizer")
    pst.parameter_data.loc[:, "partrans"] = "fixed"
    pst.parameter_data.loc[use_pars, "partrans"] = "log"
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_base.pst", num_slaves=11,
                               master_dir=test_d + "_base", verbose=True, slave_root=model_d,
                               port=4020)
    phi_df2 = pd.read_csv(os.path.join(test_d + "_base", "pest_base.phi.composite.csv"))
    par_df2 = pd.read_csv(os.path.join(test_d + "_base", "pest_base.{0}.par.csv".format(pst.control_data.noptmax)),
                          index_col=0)
    par_df2.columns = par_df1.columns.str.lower()
    plt.plot(phi_df1.total_runs, phi_df1.loc[:, "mean"], label="local")
    plt.plot(phi_df2.total_runs, phi_df2.loc[:, "mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d, "local_test.pdf"))

    for par in par_df1.columns:
        diff = par_df1.loc[:, par] - par_df2.loc[:, par]
        print(diff)
        print(diff.sum())
        diff = par_df1.loc[:, par] - pe.loc[:, par]
        print(diff.sum())
        #assert diff.sum() == 0.0
    diff = phi_df1.loc[:,"mean"] - phi_df2.loc[:,"mean"]
    print(diff.max().max())
    assert np.abs(diff.max().max()) < 0.1
  

def freyberg_localizer_test1():

    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "test_local1")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    par = pst.parameter_data
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat.loc[:,:] = 1.0
    future_pars = par.loc[par.pargp.apply(lambda x: x in ["w1","r1"]),"parnme"]
    mat.loc[:,future_pars] = 0.0
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d,"par.csv"))
   
    
    pst.observation_data.loc[pst.nnz_obs_names[3],"obgnme"] = "less_than"
    pst.observation_data.loc[pst.nnz_obs_names[3],"weight"] = 100.0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 3
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 3
    print("writing pst")
    pst.write(os.path.join(template_d, "pest_base.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_base.pst", num_slaves=11, master_dir=test_d,
                               slave_root=model_d,port=4020)
    par_df = pd.read_csv(os.path.join(test_d,"pest_base.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    par_df.index = pe.index
    par_df.columns = par_df.columns.str.lower()

    par_df_org = pd.read_csv(os.path.join(test_d, "pest_base.0.par.csv"), index_col=0)
    par_df_org.index = pe.index
    par_df_org.columns = par_df_org.columns.str.lower()

    broke = []
    for pg in ["r1","w1"]:
        o = par_df_org.loc[:,par.loc[par.pargp==pg,"parnme"]]

        omn,ostd = o.mean(),o.std()
        u = par_df.loc[:,par.loc[par.pargp==pg,"parnme"]]
        for col in o.columns:
            diff = o.loc[:,col] - u.loc[:,col]
            ad = np.abs(diff).sum()
            if ad > 1.0e-7:
                print(col,ad)
                broke.append(col)
        umn,ustd = u.mean(),u.std()
    if len(broke) > 0:
        raise Exception("future pars too diff:{0}".format(','.join(broke)))


def freyberg_localizer_test2():
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "test_local2")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    par = pst.parameter_data
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names, pst.adj_par_names).to_dataframe()
    mat.loc[:, :] = 1.0
    future_pars = par.loc[par.pargp.apply(lambda x: x in ["w1", "r1"]), "parnme"]
    mat.loc[:, future_pars] = 0.0
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par.csv"))

    pst.observation_data.loc[pst.nnz_obs_names[3], "obgnme"] = "less_than"
    pst.observation_data.loc[pst.nnz_obs_names[3], "weight"] = 100.0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 3
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 6
    print("writing pst")
    pst.write(os.path.join(template_d, "pest_local.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11, master_dir=test_d,
                               slave_root=model_d, port=4020)
    par_df1 = pd.read_csv(os.path.join(test_d, "pest_local.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df1.index = pe.index
    par_df1.columns = par_df1.columns.str.lower()
    phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))

    pst.pestpp_options.pop("ies_localizer")
    pst.write(os.path.join(template_d, "pest_base.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_base.pst", num_slaves=11, master_dir=test_d+"_base",
                               slave_root=model_d, port=4020)
    par_df2 = pd.read_csv(os.path.join(test_d+"_base", "pest_base.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df2.index = pe.index
    par_df2.columns = par_df2.columns.str.lower()
    phi_df2 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))
    plt.plot(phi_df1.total_runs, phi_df1.loc[:, "mean"], label="local")
    plt.plot(phi_df2.total_runs, phi_df2.loc[:, "mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d, "local_test.pdf"))



def freyberg_localizer_test3():
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "test_local3")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    par = pst.parameter_data
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names, pst.par_groups).to_dataframe()
    mat.loc[:, :] = 1.0
    #future_pars = par.loc[par.pargp.apply(lambda x: x in ["w1", "r1"]), "parnme"]
    mat.loc[:, ["w1","r1"]] = 0.0
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par.csv"))

    pst.observation_data.loc[pst.nnz_obs_names[3], "obgnme"] = "less_than"
    pst.observation_data.loc[pst.nnz_obs_names[3], "weight"] = 100.0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 3
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.control_data.noptmax = 6
    print("writing pst")
    pst.write(os.path.join(template_d, "pest_local.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_local.pst", num_slaves=11, master_dir=test_d,
                               slave_root=model_d, port=4020)
    par_df1 = pd.read_csv(os.path.join(test_d, "pest_local.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df1.index = pe.index
    par_df1.columns = par_df1.columns.str.lower()
    phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))

    pst.pestpp_options.pop("ies_localizer")
    pst.write(os.path.join(template_d, "pest_base.pst"))
    print("starting slaves")
    pyemu.helpers.start_slaves(template_d, exe_path, "pest_base.pst", num_slaves=11, master_dir=test_d+"_base",
                               slave_root=model_d, port=4020)
    par_df2 = pd.read_csv(os.path.join(test_d+"_base", "pest_base.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df2.index = pe.index
    par_df2.columns = par_df2.columns.str.lower()
    phi_df2 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))
    plt.plot(phi_df1.total_runs, phi_df1.loc[:, "mean"], label="local")
    plt.plot(phi_df2.total_runs, phi_df2.loc[:, "mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d, "local_test.pdf"))

def compare_freyberg_local3():
    df_l = pd.read_csv(os.path.join("ies_freyberg","test_local3","pest_local.6.par.csv"))
    df_f = pd.read_csv(os.path.join("ies_freyberg","test_local3_base","pest_base.6.par.csv"))

    pst = pyemu.Pst(os.path.join("ies_freyberg","template","pest.pst"))
    df_l.columns = df_l.columns.str.lower()
    df_f.columns = df_f.columns.str.lower()
    groups = ["w1","r1"]
    par = pst.parameter_data
    # pars = par.loc[par.pargp.apply(lambda x: x in groups)]
    # df_l = df_l.loc[:,pars]
    # df_f = df_f.loc[:,pars]
    for group,title in zip(groups,["future pumping","future recharge"]):
        fig = plt.figure(figsize=(10,10))
        ax1,ax2 = plt.subplot(211), plt.subplot(212)
        fig.suptitle(title+ " parameters with (blue) and w/o (green) localization")
        ax1.set_title("mean")
        ax2.set_title("standard deviation")

        for ax in zip(groups,[ax1,ax2]):
            pars = par.loc[par.pargp==group,"parnme"]
            df_lg = df_l.loc[:, pars]
            df_fg = df_f.loc[:, pars]
            mn_l, st_l = df_lg.mean(),df_lg.std()
            mn_f, st_f = df_fg.mean(), df_fg.std()

            print(mn_l,mn_f)
            print(st_l,st_f)
            mn_l.hist(ax=ax1,facecolor='b',alpha=0.5,normed=True)
            mn_f.hist(ax=ax1,facecolor='g',alpha=0.5,normed=True)
            st_l.hist(ax=ax2, facecolor='b', alpha=0.5,normed=True)
            st_f.hist(ax=ax2, facecolor='g', alpha=0.5,normed=True)
        ax2.set_yticklabels([])
        ax1.set_yticklabels([])

        plt.show()


if __name__ == "__main__":
    # write_empty_test_matrix()
    #setup_suite_dir("ies_freyberg")
    # setup_suite_dir("ies_10par_xsec")
    # run_suite("ies_freyberg")
    # run_suite("ies_10par_xsec")
    # rebase("ies_freyberg")
    # rebase("ies_10par_xsec")
    # compare_suite("ies_10par_xsec")
    # compare_suite("ies_freyberg")

    #tenpar_subset_test()
    #tenpar_full_cov_test()
    #test_freyberg_full_cov_reorder()
    #test_freyberg_full_cov_reorder_run()
    #test_freyberg_full_cov()
    #tenpar_tight_tol_test()
    #test_synth()
    #test_10par_xsec()
    #test_freyberg()
    #test_chenoliver()
    #tenpar_weight_pareto_test()
    #compare_pyemu()
    #tenpar_narrow_range_test()
    #test_freyberg_ineq()
    #tenpar_localizer_test1()
    #tenpar_localizer_test3()
    #freyberg_localizer_test2()
    freyberg_localizer_test3()
    #compare_freyberg_local3()
    # # invest()
    #compare_suite("ies_10par_xsec")
    #compare_suite("ies_freyberg")
    
    #test_kirishima()

    #tenpar_fixed_test()

    #setup_rosenbrock()