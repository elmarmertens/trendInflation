#!/bin/bash

source /opt/intel/oneapi/setvars.sh

export OMP_NUM_THREADS=4

caffeinate -iw $$ &

this=mcmcPaddingtonGAPSVeqf
Tdata=0

for datalabel in INF INFTRM INFTRMSRV
do
    echo $datalabel
    make THIS=$this run DATALABEL=$datalabel
done
