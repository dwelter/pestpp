import os
import shutil
import time
import numpy as np
import pandas
import pst_handler as phand

class parameter:
    def __init__(self,name,lnum,start,end,value,marker='~'):
        if len(name) > (end-start-2):
            raise Exception('not enough space to write name between markers')
        self.name = name
        self.lnum = lnum
        self.start = start
        self.end = end
        self.marker = marker
        self.value = value
    
    def insert_line(self,line):
        str_line = line[:self.start]
        end_line = line[self.end:]
        marker_string = self.build_marker_string()
        new_line = str_line + marker_string + end_line
        return new_line

    def build_marker_string(self):
        string = self.name
        while len(string) < (self.end-self.start-2):
               string = string + ' '        
        string = self.marker + string + self.marker
        return string



def write_template_file(template_filename,parameters,par_per_line=10):    
    f_tpl = open(template_filename,'w',0)    
    f_tpl.write('PTF '+parameters[0].marker+'\n')
    line = ''
    for ip,p in enumerate(parameters):
        if (ip+1) % par_per_line == 0:
            #print len(line)
            f_tpl.write(line+'\n')        
            line = ''
        line += ' ' + p.build_marker_string()
        #line = line + '\n'
        #line = p.insert_line(line)
    f_tpl.write(line+'\n')
    f_tpl.close()     


def write_par_file(par_filename,pt,precs,parameters):
    f = open(par_filename,'w',0)
    f.write(precs+' '+pt+'\n')
    for p in parameters:
        f.write(p.name+ ' {0:23.17E}  {1:1.8F}  {2:1.8F}\n'.format(p.value,1.0,0.0))
    f.close()

def write_insfile(insfile,outfile,obs_prefix,nobs):
    rand_vals = np.random.rand(nobs)
    width = 15
    starting_cols = np.random.randint(20,100,size=nobs)
    fi = open(insfile,'w')
    fo = open(outfile,'w')
    fi.write("pif $\n")
    onames = []
    for iobs in xrange(nobs):
        scol = starting_cols[iobs]
        oname = obs_prefix+"{0:d}".format(iobs)
        onames.append(oname)
        name_fmt = "{0:>"+str(scol)+"s}"
        fo.write(name_fmt.format(oname)+"{0:15.6E}\n".format(rand_vals[iobs]))
        fi.write("l1 ["+oname+"]{0:d}:{1:d}\n".format(scol+1,scol+1+width))
    fo.close()
    fi.close()
    return onames


exe_path = '..\\x64\\Release\\'
#exe_path = '..\\Debug\\'
exe_name = 'pest++.exe'

val = 1.12345678901234567890e15
        
scale = 1.0
offset = 0.0
f_log = open('iopp.log','w',0)
f_log.write('running tests...\n')

test_count = 10
tdir = 'test_results\\'
npar = 20000
nobs = 100000

shutil.rmtree(tdir,ignore_errors=True)
try:
    os.mkdir(tdir)
except:
    os.mkdir(tdir)

test_failed = []


start = 0
end = 10
parameters = []
pnames = []                                                   
for pnum in range(npar):
    pname = 'p'+str(pnum+1)
    pnames.append(pname)
    par = parameter(pname,pnum+1,start,end,val)
    parameters.append(par)      
ins_file = os.path.join(tdir,"test.ins")  
out_file = os.path.join(tdir,"test.bak.out")     
tpl_file = os.path.join(tdir,"test.tpl")                    
write_template_file(tpl_file,parameters)
obs_names = write_insfile(ins_file,out_file,"obs_",nobs)

pst = phand.pst("test.bak.pst")

#"PARNME PARTRANS PARCHGLIM PARVAL1 PARLBND PARUBND PARGP SCALE OFFSET DERCOM"
pvals = np.ones(len(pnames))
partrans = ["log"] * len(pnames)
parchglim = ["factor"] * len(pnames)
parlbnd,parubnd = np.zeros(len(pnames))+0.1,np.ones(len(pnames))+2.0
pargp = ["par"] * len(pnames)
scale = np.ones(len(pnames))
offest = np.zeros(len(pnames))
dercom = np.ones(len(pnames))
par = pandas.DataFrame({"parnme":pnames,"partrans":partrans,"parchglim":parchglim,\
	"parval1":pvals,"parlbnd":parlbnd,"parubnd":parubnd,"pargp":pargp,"scale":scale,\
	"offset":offset,"dercom":dercom},index=pnames,columns=pst.par_fieldnames)
pst.parameter_data = par
print pst.parameter_data.columns
#"OBSNME OBSVAL WEIGHT OBGNME"
weight = np.ones(len(obs_names))
obsval = np.ones(len(obs_names))
obgnme = ["obs"] * len(obs_names)
obs = pandas.DataFrame({"obsnme":obs_names,"obsval":obsval,"weight":weight,"obgnme":obgnme},index=obs_names,columns=pst.obs_fieldnames)
pst.observation_data = obs
pst.write(os.path.join(tdir,"test.pst"))
f = open(os.path.join(tdir,"model.bat"),'w')
f.write("copy test.bak.out test.out\n")
f.close()


for itest in xrange(test_count):
	

