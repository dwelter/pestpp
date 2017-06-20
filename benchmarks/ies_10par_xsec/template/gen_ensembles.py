import pyemu

pst = pyemu.Pst("pest.pst")
mc = pyemu.MonteCarlo(pst=pst)
mc.draw(5,obs=True)
mc.parensemble.to_csv("par.csv")
mc.obsensemble.to_csv("obs.csv")