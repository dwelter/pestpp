import pyemu
pst = pyemu.Pst("dewater_pest.pst")
jco = pyemu.Jco.from_binary("dewater_pest.jcb")
print(jco.x)
oc_dict = {oc:"l" for oc in pst.nnz_obs_names}
obj_func = {name:1.0 for name in pst.par_names}
pyemu.utils.to_mps("dewater_pest.jcb",obj_func=obj_func,obs_constraint_sense=oc_dict,pst=pst)



