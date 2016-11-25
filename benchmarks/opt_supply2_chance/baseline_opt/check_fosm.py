import numpy as np
import pyemu

pst = pyemu.Pst("supply2_pest.pst")
base_par_names = pst.adj_par_names

new_par_names = [name for name in base_par_names if name[0] in ["t","s"]]



print(new_par_names)

new_pst = pst.get(par_names=new_par_names)

forecast_names = pst.observation_data.groupby(pst.observation_data.obgnme).groups["less_obs"]

sc = pyemu.Schur(jco="supply2_pest.3.jcb",pst=new_pst,forecasts=forecast_names)

sum = sc.get_forecast_summary().loc[:,"post_var"].apply(np.sqrt)
print(sum)