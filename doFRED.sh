#!/bin/bash

initifort

export OMP_NUM_THREADS=4

this=mcmcPaddingtonGAPSVeqf
Tdata=0

for datalabel in INF INFTRM INFTRMSRV
do
    echo $datalabel
    make THIS=$this run DATALABEL=$datalabel
done
