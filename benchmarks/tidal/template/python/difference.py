import os
from datetime import datetime
import numpy as np
import pandas
'''makes huge assumptions about the sorting of the dataframes
'''
base_group = "O10"

def parse_date(x):
    return datetime.strptime(x,"%d/%m/%Y %H:%M:%S")

def time_series():
    in_smp = os.path.join("smp","mod","hds.smp")
    out_smp = os.path.join("smp","mod","hds_diff.smp")    
    f_out = open(out_smp,'w')
    df = pandas.read_csv(in_smp,date_parser=parse_date,\
            parse_dates=[[1,2]],sep="\\s+",header=None,\
            names=["site","date","time","value"],\
            na_values=["dry_or_inactive"])            
    site_groups = df.groupby("site").groups
    base_df = df.ix[site_groups[base_group]]
    
    for site,idxs in site_groups.iteritems():
        if site != base_group:
            name = site+"_td"
            site_df = df.ix[idxs]
            diff = base_df["value"].values - site_df["value"].values            
            for dt,val in zip(base_df["date_time"],diff):
                f_out.write("{0:20s} ".format(name)+"  "+dt.strftime("%d/%m/%Y %H:%M:%S")+" {0:15.6E}\n".format(val))
    f_out.close()                
    
if __name__ == "__main__":    
    time_series()                
