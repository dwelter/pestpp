#!/bin/bash
rm -f model/ref/hk_layer_1.ref
exe/fac2real < misc/fac2real.in >> /dev/null

cd model
./mf2005 syn.nam >> /dev/null
cd ..

exe/mod2smp < misc/mod2smp_hds.in >> /dev/null
