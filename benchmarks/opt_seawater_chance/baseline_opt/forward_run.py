import os
import numpy as np
hds_file = 'seawater.hds'
nrow,ncol,nlay = 20,30,2
delt = 200
grad_file = 'gradients.dat'
def process_hds():

    arr = np.loadtxt(hds_file).reshape((nlay,nrow,ncol))
    hg1 = (arr[0,:,27] - arr[0,:,29]) / (2.0 * delt)
    hg2 = (arr[1,:,27] - arr[0,:,29]) / (2.0 * delt)
    vg1 = (arr[0,:,28] - arr[1,:,28])
    print(arr[0,:,28], arr[1,:,28])
    hg1_labels = ["hg1_r{0}".format(i+1) for i in range(nrow)]
    hg2_labels = ["hg2_r{0}".format(i+1) for i in range(nrow)]
    vg1_labels = ["vg1_r{0}".format(i+1) for i in range(nrow)]
    with open(grad_file,'w') as f:
        f.write("obsnme obsval\n")
        for arr,labels in zip([hg1,hg2,vg1],[hg1_labels,hg2_labels,vg1_labels]):
            for v,l in zip(arr,labels):
                f.write("{0} {1:15.6E}\n".format(l,v))


if __name__ == '__main__':
    os.system('mf2005 seawater.nam')
    process_hds()
