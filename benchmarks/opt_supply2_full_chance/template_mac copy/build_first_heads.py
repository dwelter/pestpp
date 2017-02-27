import pyemu

pst = pyemu.Pst("supply2_pest.pst")
obs = pst.observation_data
obs.loc[:,"weight"] = 0.0
first_names = obs.groupby(obs.obsnme.apply(lambda x: x.endswith("_01"))).groups[True]
print(first_names[::4])
obs.loc[first_names[::4],"weight"] = 1.0
par = pst.parameter_data
t_unfixed = [name for name in par.parnme[::20] if name.startswith('t')]
par.loc[t_unfixed,"partrans"] = "log"
pst.pestpp_options["base_jacobian"] = "supply2_pest.full.jcb"
pst.pestpp_options["opt_recalc_fosm_every"] = 100
pst.write("supply2_pest.first.pst")