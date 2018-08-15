#!/bin/bash

DTF_DIR= #set path to where the DTF library is
PNETCDF_DIR= #set path to the installed PnetCDF library

export LD_LIBRARY_PATH=$DTF_DIR:$PNETCDF_DIR/lib:$LD_LIBRARY_PATH
export DTF_VERBOSE_LEVEL=0
export MAX_WORKGROUP_SIZE=2
export DTF_GLOBAL_PATH=.

px=40
py=40
pz=20

psx=2
psy=2
psz=1

nproc=$(($psx*$psy*$psz))
mpiexec -n $nproc ./writer/writer $px $py $pz $psx $psy $psz 1 F . &
mpiexec -n $nproc ./reader/reader $px $py $pz $psx $psy $psz 1 T .
wait

