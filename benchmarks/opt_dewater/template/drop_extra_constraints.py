import pyemu

pst = pyemu.Pst("dewater_pest.pst")
nnz_obs = pst.nnz_obs_names
print(nnz_obs)
ins_name = "dewater_pest.hds.drop.ins"
pst = pst.get(obs_names=nnz_obs)
pst.instruction_files = [ins_name]
pst.write("dewater_pest.drop.pst")

nrow,ncol = 20,30

with open(ins_name,'w') as f:
    f.write("pif ~\n")
    for i in range(nrow):
        if i == 0:
            f.write("l2")
        else:
            f.write("l1")
        for j in range(ncol):
            name = "h{0:02d}_{1:02d}".format(i,j)
            if name in nnz_obs:
                name = '!{0}!'.format(name)
            else:
                name = "!dum!"
            f.write(' w {0}'.format(name))
        f.write('\n')