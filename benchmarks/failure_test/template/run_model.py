import time
from subprocess import call

# read input file
fin = open("input.dat")
line = fin.next()
atoms= line.strip().split()
rech = float(atoms[0])
cond = float(atoms[1])
scoeff = float(atoms[2])
if scoeff == .15:
    print 'sleeping'
    time.sleep(60*30)  #sleep for 30 min
elif scoeff == .18:
    print 'failed run'
    pass
else:
    print 'running model'
    call("storage.exe")
