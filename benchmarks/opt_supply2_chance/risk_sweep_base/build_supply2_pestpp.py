import os
import numpy as np
import pandas as pd
import flopy
import pyemu

nrow = 25
ncol = 30
    

def build_wel_file():
    #wel_data = [[1, 12, 11, 0.0, "Q1"], [1, 16, 17, 0.0, "Q2"], [1, 11, 22, 0.0, "Q3"], [1, 14, 25, 0.0, "Q4"]]
    wel_data = [[1, 12, 11, 0.0, "Q1"], [1, 16, 17, 0.0, "Q2"], [1, 14, 25, 0.0, "Q4"]]
    ml = flopy.modflow.Modflow.load("supply2.nam", check=False, verbose=True)
    wel_data = {i: wel_data for i in range(ml.nper)}
    dtype = np.dtype([("layer", int), ("row", int), ("column", int),
                      ("flux", np.float32), ("name", object)])
    wel = flopy.modflow.ModflowWel(ml, stress_period_data=wel_data, dtype=dtype)
    wel.write_file()

def write_command_script():
    with open("supply2.py",'w') as f:
        f.write("import os\n")
        f.write("import depletion_constraints as dc\n")
        f.write("os.system('mf2005.exe supply2.nam')\n")
        f.write("dc.apply()\n")

def write_external_par_tpl():
    with open("external.tpl",'w') as f:
        f.write("ptf ~\n")
        for pname in ["im9","im10","im11","im12"]:
            f.write(" ~   " + pname + "      ~\n")

def write_tran_tpl():
    with open("TRAN.ref",'w') as f_in, open("TRAN.ref.tpl",'w') as f_tpl:
        f_tpl.write("ptf ~\n")
        for i in range(nrow):
            for j in range(ncol):
                pname = "t{0:02d}_{1:02d}".format(i,j)
                f_tpl.write(" ~   {0}  ~".format(pname))
            f_tpl.write("\n")
    np.savetxt("TRAN.ref",np.zeros((nrow,ncol))+500.0,fmt="%15.6E",delimiter='')

def write_hds_instruction_file():
    with open('supply2.hds.ins','w') as f:
        f.write("pif ~\n")
        for i in range(nrow):
            if i == 0:
                f.write("l2 ")
            else:
                f.write("l1 ")
            for j in range(ncol):
                name = "!h{0:02d}_{1:02d}!".format(i,j)
                f.write(" w {0}".format(name))
            f.write("\n")



def build_control_file():
    
    tpl_files = ["supply2.wel.tpl",'external.tpl',"supply2.sfr.tpl","tran.ref.tpl"]
    in_files = ["supply2.wel","external.dat","supply2.sfr","tran.ref"]
    ins_file = ["sfr_aq_ex.dat.ins","supply2.hds.ins"]
    out_file = ["sfr_aq_ex.dat","supply2.hds"]
    pst = pyemu.pst_utils.pst_from_io_files(tpl_files,in_files,ins_file,out_file)
    par = pst.parameter_data
    par.loc[:,"partrans"] = "none"
    par.loc[:,"parval1"] = 0.0
    par.loc[:,"parlbnd"] = 0.0
    par.loc[:,"parubnd"] = 50000.0
    par.loc[:,"scale"] = -1.0
    par.loc[:,"pargp"] = "pumping"
    im_names = par.groupby(par.parnme.apply(lambda x:x.startswith("im"))).groups[True]
    par.loc[im_names,"pargp"] = "external"
    par.loc[im_names,"parubnd"] = 1.0E+6

    t_names = par.groupby(par.parnme.apply(lambda x: x.startswith("t"))).groups[True]
    par.loc[t_names, "pargp"] = "trans"
    par.loc[t_names,"parval1"] = 5000.0
    par.loc[t_names, "parubnd"] = 6000.0
    par.loc[t_names, 'parlbnd'] = 4000.0
    par.loc[t_names,"partrans"] = "log"
    par.loc[t_names, "scale"] = 1.0

    seg_names = par.groupby(par.parnme.apply(lambda x: x.startswith("seg"))).groups[True]
    par.loc[seg_names, "pargp"] = "sfr_seg"
    par.loc[seg_names, "parval1"] = 5.0
    par.loc[seg_names, "parubnd"] = 5.5
    par.loc[seg_names, 'parlbnd'] = 4.5
    par.loc[seg_names, "partrans"] = "log"
    par.loc[seg_names, "scale"] = 1.0

    pst._rectify_pgroups()

    pst.parameter_groups.loc["pumping","inctyp"] = "absolute"
    pst.parameter_groups.loc["pumping","derinc"] = 5000.0

    const_dict = {}
    const_dict["s1r14_09"] = 15000.0
    const_dict["s1r14_10"] = 15000.0
    const_dict["s1r14_11"] = 15000.0
    const_dict["s1r14_12"] = 15000.0
    const_dict["s1r21_11"] = 20000.0
    const_dict["s1r21_12"] = 20000.0
    const_dict["s2r08_09"] = 15000.0
    const_dict["s2r08_10"] = 15000.0
    const_dict["s2r08_11"] = 15000.0
    const_dict["s2r08_12"] = 15000.0
    const_dict["s3r05_11"] = 30000.0
    const_dict["s3r05_12"] = 30000.0
    obs = pst.observation_data
    obs.loc[:,"obgnme"] = "extra_sfr"
    obs.loc[:,"obsval"] = 0.0
    obs.loc[:,"weight"] = 0.0
    obs.loc[const_dict.keys(),"obgnme"] = "less_obs"
    for obsnme,value in const_dict.items():
        obs.loc[obsnme,"obsval"] = value
        obs.loc[obsnme,"weight"] = 1.0

    h_groups = obs.groupby(obs.obsnme.apply(lambda x:x.startswith('h'))).groups[True]
    obs.loc[h_groups,"obgnme"] = "head"

    obs.sort(["obgnme","obsnme"],inplace=True)

    pst.control_data.noptmax = 1
    pst.model_command = ["python.exe supply2.py"]
    obj_eq = ' 1.095 * q1 + 0.276 * q2a + 0.273 * q2b + 0.273 * q2c + 0.273 * q2d + 0.273 * q4a + 0.273 * q4b'
    obj_eq += ' + -0.1104 + im9 + -0.1092 * im10 + -0.01092 * im11 + -0.1092 * im12 = 0.0'
    # for pname in pst.adj_par_names:
    #     if pname.startswith("im"):
    #         obj_eq += " + -0.0012 * " + pname
    #         #obj_eq += " + 0.0012 * " + pname
    #     elif pname.startswith("q"):
    #         obj_eq += " + 0.001 * " + pname
    #         #obj_eq += " + -0.001 * " + pname
    # obj_eq += ' = 0.0'
    # obj_eq = obj_eq[2:]
    #pst.prior_information.loc["obj_func","pilbl"] = "obj_func"
    #pst.prior_information.loc["obj_func","equation"] = obj_eq
    #pst.prior_information.loc["obj_func","weight"] = 1.0
    #pst.prior_information.loc["obj_func","obgnme"] = "e"

    pi_eqs,pi_weights,pi_names,pi_groups = [obj_eq],[],["obj_func"],['obj_func']
    #sp0
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a = 30000.0")
    pi_groups.append("greater_pi")
    
    #sp1
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a = 30000.0")
    pi_groups.append("greater_pi")
    
    #sp2
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c = 30000.0")
    pi_groups.append("greater_pi")
    
    #sp3
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b = 30000.0")
    pi_groups.append("greater_pi")
    
    #sp4
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a = 25000.0")
    pi_groups.append("greater_pi")
    
    #sp5
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a = 25000.0")
    pi_groups.append("greater_pi")
    
    #sp6
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c = 25000.0")
    pi_groups.append("greater_pi")
    
    #sp7
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b = 25000.0")
    pi_groups.append("greater_pi")

    #sp8
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a + 1.0 * im9 = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2a + 1.0 * im9 = 45000.0")
    pi_groups.append("greater_pi")
    
    #sp9
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a + 1.0 * im10 = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2b + 1.0 * q4a + 1.0 * im10  = 45000.0")
    pi_groups.append("greater_pi")
    
    #sp10
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c + 1.0 * im11 = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2c + 1.0 * im11 = 45000.0")
    pi_groups.append("greater_pi")
    
    #sp11
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b + 1.0 * im12 = 80000.0")
    pi_groups.append("less_pi")
    pi_eqs.append(" 1.0 * q1 + 1.0 * q2d + 1.0 * q4b + 1.0 * im12 = 45000.0")
    pi_groups.append("greater_pi")

    pi_weights = [1.0 for _ in range(len(pi_eqs))]
    pi_names = ["pi_const{0}".format(i) for i in range(len(pi_eqs))]
    pi_names[0] = "pi_obj_func"

    #print(len(pi_eqs),len(pi_weights),len(pi_names),len(pi_groups))

    pst.prior_information = pd.DataFrame({"pilbl":pi_names,"equation":pi_eqs,
        "weight":pi_weights,"obgnme":pi_groups},index=pi_names)

    pst.pestpp_options["opt_direction"] = "max"
    pst.pestpp_options["opt_constraint_groups"] = "less_pi,greater_pi,less_obs"
    pst.pestpp_options["opt_obj_func"] = "pi_obj_func"
    pst.pestpp_options["opt_coin_log"] = "true"
    pst.pestpp_options["opt_dec_var_groups"] = "pumping,external"
    pst.pestpp_options["opt_ext_var_groups"] = "external"
    #pst.pestpp_options["base_jacobian"] = "supply2_pest.1.bak.jcb"
    pst.parameter_data.sort_index(inplace=True)
    pst.write("supply2_pest.pst")
    print(pst.nnz_obs_names)


def get_obj_func():
    dv_counts = {}
    with open("supply2.wel.tpl",'r') as f:
        f.readline()
        while True:
            line = f.readline()
            if line == '':
                break
            if "stress period" in line.lower():
                #print(line)
                nwel = int(line.split()[0])
                iper = int(line.split()[-1])
                for _ in range(nwel):
                    line = f.readline()
                    if line == '':
                        raise exception()
                    if "~" in line:
                        w = line.split('~')[1].split()[0]
                        if w not in dv_counts:
                            dv_counts[w] = []
                        dv_counts[w].append(iper)
    ml = flopy.modflow.Modflow.load("supply2.nam",load_only=[])
    perlen = ml.dis.perlen.array
    dv_coefs = {}
    rate = 0.001

    for dv,pers in dv_counts.items():
        dv_coefs[dv] = "{0:5.4f}".format(rate * perlen[pers].sum())
    #print(dv_coefs)

    im_coefs = {}
    im_rate = 0.0012
    for iper in range(9,13):
        im_coefs["im{0}".format(iper)] = "{0:5.4f}".format(im_rate * perlen[iper-1])
    
    obj_str = ''
    dvs = list(dv_coefs.keys())
    dvs.sort()
    for dv in dvs:
        obj_str += ' {0} * {1} +'.format(dv,dv_coefs[dv])
    obj_str = obj_str[:-1]
    #print(obj_str)
    
    ims = list(im_coefs.keys())
    ims.sort()
    for im in ims:
        obj_str += ' - {0} * {1}'.format(im,im_coefs[im])
    obs_str = obj_str[:-1] + ' = 0.0'
    #print(obj_str)
    return obj_str


if __name__ == "__main__":
    # build_wel_file()
    #write_external_par_tpl()
    #write_command_script()
    #write_tran_tpl()
    #write_hds_instruction_file()
    #build_control_file()
    get_obj_func()