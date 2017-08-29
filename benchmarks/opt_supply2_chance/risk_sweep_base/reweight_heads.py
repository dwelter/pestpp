import pyemu

pst = pyemu.Pst("supply2_pest.full.pst")
obs = pst.observation_data

head_names = obs.groupby("obgnme").groups["head"]

obs.loc[head_names[::3],"weight"] = 1.0

sfr_1_names = obs.groupby(obs.apply(lambda x: x.obsnme.startswith("s1") and x.obgnme.startswith("extra"),axis=1)).groups[True]

obs.loc[sfr_1_names[::3],"weight"] = 1.0

pst.write("supply2_pest.some.pst")