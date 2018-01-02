import os
import pyemu

gs = pyemu.utils.read_struct_file(os.path.join("misc","struct.dat"))[0]
pp_df = pyemu.gw_utils.pp_file_to_dataframe(os.path.join("misc","pp_locs.dat"))
cov = gs.covariance_matrix(x=pp_df.x,y=pp_df.y,names=pp_df.name)
cov.to_ascii(os.path.join("misc","pp.cov"))