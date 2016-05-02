import os
import shutil
from datetime import datetime
import pst_handler as phand
import tsprocClass as tc
import pestUtil as pu

tc.DATE_FMT = "%d/%m/%Y"

#start,end = datetime(year=1999,month=1,day=1),datetime(year=2009,month=12,day=31)

def main():
    
    #--instance
    tsproc_infile = 'tsproc_setup.dat'
    tsp = tc.tsproc(tsproc_infile,out_file='processed.dat',out_fmt='long',)
    if os.path.exists('pest.pst'):
        os.remove('pest.pst')
    pest_mblocks,pest_oblocks = [],[]

   
    

    #---------------------------
    #--heads
    #--------------------------
    print 'procesing heads'
    ohead_name = 'smp\\obs\\hds.smp'
    mhead_name = 'smp\\mod\\hds.smp'
    #shutil.copy2(mhead_name,ohead_name)
    
    
    ohead = pu.smp(ohead_name,load=False,site_index=True)    
    mhead = pu.smp(mhead_name,load=False,site_index=True)
    #--some model sites are not included - thin them out

    keep_mlist = []
    names = []
    mlist = list(mhead.site_list)
    olist = list(ohead.site_list)
    for site in mlist:
        if site in olist:
            keep_mlist.append(site)
        else:
            print 'missing mhead site',site
            #missing.append(site)


    oblocks = tsp.get_mul_series_ssf(keep_mlist,ohead_name,suffix='oh',\
        role="final",context=tc.PEST_CONTEXT)
    mblocks = tsp.get_mul_series_ssf(keep_mlist,mhead_name,suffix='sh')
    renamed = tsp.copy_2_series(mblocks,keep_mlist,role='final',wght=1.0)
    pest_oblocks.extend(oblocks)
    pest_mblocks.extend(renamed)

   
    print 'build a list of template and model-equivalent files'
    tpl_dir = 'tpl\\'
    modin_dir = os.path.join("model","ref")
    tpl_files,modin_files = [],[]
    files = os.listdir(tpl_dir)
    for f in files:
        modin_file = f.replace(".tpl",'')
        modin_files.append(os.path.join(modin_dir,modin_file))
        tpl_file = os.path.join(tpl_dir,f)        
        tpl_files.append(tpl_file)
            
    par_file = os.path.join("misc","par.info")
    grp_file = os.path.join("misc","grp.info")

    #--write the model run tspoc file 
    print 'writing model run tsproc file'                             
    tsp.set_context('model_run')
    tsp.tsproc_file = 'misc\\tsproc_model_run.dat'
    tsp.write_tsproc()

    #--write the setup tsproc file
    print 'writing setup tsproc file'
    tsp.write_pest(tpl_files,modin_files,pest_oblocks,pest_mblocks,\
                   svd=True,parms=par_file,parm_grp=grp_file,\
                   model_cmd="model.bat")
    
    tsp.set_context(tc.PEST_CONTEXT)
    tsp.tsproc_file = 'misc\\tsproc_setup.dat'
    tsp.write_tsproc()

    f = open(os.path.join('misc','tsproc_setup.in'),'w')
    f.write(os.path.join('misc','tsproc_setup.dat')+'\n'+\
                        os.path.join('misc','tsproc_setup.out')+\
                            '\ny\ny\n')
    f.close()

    f = open(os.path.join('misc','tsproc_model_run.in'),'w')
    f.write(os.path.join('misc','tsproc_model_run.dat')+'\n'+\
            os.path.join('misc','tsproc_model_run.out')+'\ny\ny\n')
    f.close()

    print 'running tsproc'
    os.system(os.path.join("exe",'tsproc.exe')+' <'+os.path.join("misc","tsproc_setup.in")+\
              ' >'+os.path.join("misc","tsproc_screen.out"))


if __name__ == '__main__':
    main()



