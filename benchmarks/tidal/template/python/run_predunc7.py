import os
import pandas 
import pst_handler as phand

pst_dir = os.path.join("pst_jco")
base_case = "pest"

files = os.listdir(pst_dir)
psts = []
for f in files:
    if f.endswith(".pst") and f.split('.')[0] != base_case:
        psts.append(f)
print psts  

for pst in psts:
    in_case = os.path.join(pst_dir,base_case)
    out_case = os.path.join(pst_dir,pst.split('.')[0])
    exe = os.path.join("exe","i64jco2jco.exe")
    os.system(" ".join([exe,in_case,out_case]))
    
    args = [os.path.join(pst_dir,pst),"1.0",os.path.join("unc","unc.dat"),\
            os.path.join("unc",pst+".cov"),"1"]    
    in_name = os.path.join("unc",pst+".in")             
    f = open(in_name,'w')
    f.write('\n'.join(args)+'\n')
    f.close()    
    exe = os.path.join("exe","predunc7.exe")
    os.system(exe + " <"+in_name)
    
 
        
