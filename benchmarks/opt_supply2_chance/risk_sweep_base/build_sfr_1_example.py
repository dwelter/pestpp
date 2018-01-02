import pyemu

pst = pyemu.Pst("supply2_pest.pst")

obs = pst.observation_data
sfr_1_names = obs.groupby(obs.obsnme.apply(lambda x:x.startswith("h") and x.endswith("_01"))).groups[True]
obs.loc[sfr_1_names,"weight"] = 1.0
pst.write("supply2_pest.sfr1.pst")