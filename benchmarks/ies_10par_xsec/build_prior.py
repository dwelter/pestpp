import os
import pandas as pd
import pyemu

#pp = pd.read_csv(os.path.join("template","misc","pp_locs.dat"),delim_whitespace=True,
#                 header=None,names=pyemu.gw_utils.PP_NAMES)
x = [float(i)*10.0 for i in range(10)]
y = [0.0] * 10
pst = pyemu.Pst(os.path.join("master_parcov","pest.pst"))
names = list(pst.parameter_data.parnme)
names.remove("stage")
ev = pyemu.utils.ExpVario(0.001,500.0)
ppcov = ev.covariance_matrix(x=x,y=y,names=names)
pstcov = pyemu.Cov.from_parbounds(os.path.join("template","pest.pst"))
pstcov.drop(names,axis=1)
ppcov = ppcov.extend(pstcov)
ppcov.to_ascii(os.path.join("master_parcov","prior.cov"))