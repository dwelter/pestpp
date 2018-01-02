import numpy as np
import pyemu

pst = pyemu.Pst("supply2_pest.sfr1.fromfulljco.pst")
obs = pst.observation_data
pst.prior_information = pst.null_prior
forecasts = obs.groupby(obs.obgnme).groups["less_obs"]
pars = [par for par in pst.adj_par_names if not par[0] in ["i",'q']]
obs = [obs for obs in pst.nnz_obs_names if obs not in forecasts]
print(pars,obs)
sc = pyemu.Schur(jco="supply2_pest.full.jcb",pst=pst,forecasts=forecasts)
sc = sc.get(par_names=pars,obs_names=obs)
print(sc.pst.npar_adj,sc.pst.nnz_obs)
df = sc.get_forecast_summary().loc[:,["prior_var","post_var"]].apply(np.sqrt)
print(df)

sc.posterior_parameter.to_ascii("supply2_pest.fosm.post.cov")
pst = pyemu.Pst("supply2_pest.sfr1.fromfulljco.pst")
obs = pst.observation_data
names = obs.groupby(obs.obgnme.apply(lambda x:x=="less_obs")).groups[False]
obs.loc[names,"weight"] = 0.0
pst.pestpp_options["parcov_filename"] = "supply2_pest.fosm.post.cov"
pst.write("supply2_pest.sfr1.fromfulljco.postcov.pst")