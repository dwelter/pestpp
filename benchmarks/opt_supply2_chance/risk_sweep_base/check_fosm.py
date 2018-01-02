import numpy as np
import pyemu

pst = pyemu.Pst("supply2_pest.fosm.pst")
base_par_names = pst.adj_par_names

new_par_names = [name for name in base_par_names if name[0] in ["t","s"]]



print(new_par_names)

new_pst = pst.get(par_names=new_par_names)

forecast_names = pst.observation_data.groupby(pst.observation_data.obgnme).groups["less_obs"]
#new_pst.observation_data.loc[forecast_names,"weight"] = 0.0
sc = pyemu.Schur(jco="supply2_pest.full.jcb",pst=new_pst,forecasts=forecast_names,verbose=True)

print(sc.obscov)

#sum = sc.get_forecast_summary().loc[:,"post_var"].apply(np.sqrt)
#print(sum)

df_worth = sc.get_added_obs_importance(reset_zero_weight=1.0)
df_worth.to_csv("worth.csv")