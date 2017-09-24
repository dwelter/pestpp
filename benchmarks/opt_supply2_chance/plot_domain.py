import os
import numpy as np
import pandas as pd
import matplotlib
font = {'size'   : 8}

matplotlib.rc('font', **font)
from matplotlib.patches import Rectangle,Polygon
import matplotlib.pyplot as plt
import flopy
import pyemu


model_ws = os.path.join("template_mac")
pst = pyemu.Pst(os.path.join(model_ws,"supply2_pest.base.pst"))
# pst.parrep(os.path.join("baseline_opt","supply2_pest.base.par"))
# pst.control_data.noptmax = 0
# pst.write(os.path.join(model_ws,"supply2_pest.base.final.pst"))
# os.chdir(model_ws)
# os.system("pestpp supply2_pest.base.final.pst")
# os.chdir("..")

obs = pst.observation_data
sfr_const = obs.loc[obs.groupby("obgnme").groups["less_obs"],:]
sfr_const.loc[:,"segment"] = sfr_const.obsnme.apply(lambda x:int(x[1]))
sfr_const.loc[:,"reach"] = sfr_const.obsnme.apply(lambda x:int(x[3:5]))
sfr_const.loc[:,"tup"] = sfr_const.apply(lambda x: (x.segment,x.reach),axis=1)
sfr_const.loc[:,"tup_str"] = sfr_const.tup.apply(lambda x: "{0}_{1}".format(*x))
const_tups = sfr_const.tup_str.unique()

ml = flopy.modflow.Modflow.load("supply2.nam",model_ws=model_ws,check=False)
hds = flopy.utils.FormattedHeadFile(os.path.join(model_ws,ml.name+".hds"),model=ml)
harray = hds.get_data()[0]

reach_data = pd.DataFrame.from_records(ml.sfr.reach_data)
reach_data.loc[:,"iseg"] = reach_data.loc[:,"iseg"].apply(np.int)
reach_data.loc[:,"ireach"] = reach_data.loc[:,"ireach"].apply(np.int)
reach_data.loc[:,"tup"] = reach_data.apply(lambda x: (int(x.iseg),int(x.ireach)),axis=1)
reach_data.loc[:,"tup_str"] = reach_data.tup.apply(lambda x:"{0}_{1}".format(*x))
reach_data.loc[:,"x"] = reach_data.apply(lambda x:ml.sr.xcentergrid[x.i,x.j],axis=1)
reach_data.loc[:,"y"] = reach_data.apply(lambda x:ml.sr.ycentergrid[x.i,x.j],axis=1)
isconst = reach_data.tup_str.apply(lambda x:x in const_tups)
reach_data_const = reach_data.loc[isconst,:]

#print(ml.sfr.reach_data.dtype)

wel_data = ml.wel.stress_period_data[0]
wel_x = [ml.sr.xcentergrid[i,j] for i,j in zip(wel_data["i"],wel_data["j"])]
wel_y = [ml.sr.ycentergrid[i,j] for i,j in zip(wel_data["i"],wel_data["j"])]
wel_label = ["q{0}".format(i+1) for i in range(len(wel_x))]

for w in [wel_x,wel_y,wel_label]:
    w.pop(2)

fig = plt.figure(figsize=(17.15/2.54,10.15/2.54))
ax = plt.axes([0.1,0.075,0.7,0.9])
cax = plt.axes([0.85,0.075,0.025,0.9])
#i = ax.imshow(harray,interpolation="nearest",alpha=0.5,extent=ml.sr.get_extent())
i = ax.pcolormesh(ml.sr.xedge,ml.sr.yedge,harray,alpha=0.5)
c = plt.colorbar(i,cax=cax)
c.set_label("water level (ft)",labelpad=0.1)

ax.set_xlabel("x distance (ft)",labelpad=0.1)
ax.set_ylabel("y distance (ft)",labelpad=0.1)
ax.scatter(wel_x,wel_y,marker='.',s=10,color='k',label="pumping well cell",zorder=10)
for x,y,l in zip(wel_x,wel_y,wel_label):
    t = ax.text(x+10,y+10,l,ha="left",va="bottom")
    t.set_bbox(dict(color='w', edgecolor='none',pad=0.0))
xmn,xmx,ymn,ymx = ml.sr.get_extent()
# for i,j in zip(ml.sfr.reach_data["i"],ml.sfr.reach_data['j']):
#     pts = np.array(ml.sr.get_vertices(i,j))
#     p = Polygon(pts,facecolor="none",edgecolor='k',lw=0.5)
#     #ax.plot(pts[:,0],pts[:,1],"k-",lw=0.5)
#     ax.add_patch(p)
# for i,j in zip(reach_data.i,reach_data.j):
#     print(i,j)
#     pts = np.array(ml.sr.get_vertices(i,j))
#     p = Polygon(pts,facecolor="0.5",alpha=0.5,edgecolor='k',lw=0.5)
#     ax.add_patch(p)
r1 = Rectangle((0,0),0,0,facecolor="none",edgecolor='b',lw=0.5,label="sfr cell")
ax.add_patch(r1)

seg_group = reach_data.groupby("iseg").groups
for iseg,idxs in seg_group.items():
    ax.plot(reach_data.loc[idxs,"x"],reach_data.loc[idxs,"y"],"k--")
ax.scatter(reach_data_const.x,reach_data_const.y,marker='^',s=15,color='k',zorder=10)
for i,(x,y) in enumerate(zip(reach_data_const.x,reach_data_const.y)):
    t = ax.text(x+10,y+10,"c{0}".format(i+1),ha="left",va="bottom")
    t.set_bbox(dict(color='w', edgecolor='none',pad=0.0))
xmn,xmx,ymn,ymx = ml.sr.get_extent()
ax.set_xlim(xmn,xmx)
ax.set_ylim(ymn,ymx)
#plt.tight_layout()
plt.savefig("supply2.pdf")




