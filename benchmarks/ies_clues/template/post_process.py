import os
import sys
import shutil
import pandas as pd

proj_default = "Project_Default"
if len(sys.argv) > 1:
	proj_default = sys.argv[1]

df = pd.read_csv(os.path.join(proj_default,"Outputs","Reaches_Output_Short.csv"),index_col=0)

df.to_csv("processed_reach_output.dat",sep=' ')