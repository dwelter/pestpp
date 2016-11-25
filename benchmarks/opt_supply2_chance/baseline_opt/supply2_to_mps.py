import pyemu
pst = pyemu.Pst("supply2_pest.pst")

oc_dict = {oc:"l" for oc in pst.nnz_obs_names}
obj_func = {name:-1.095 for name in pst.par_names}
pyemu.utils.to_mps("supply2_pest.jcb",obj_func=obj_func,obs_constraint_sense=oc_dict,pst=pst)



