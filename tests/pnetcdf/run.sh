#!/bin/bash

export DTF_GLOBAL_PATH=.
export LD_LIBRARY_PATH=../../libdtf:$HOME/apps/pnetcdf/lib:$LD_LIBRARY_PATH
export DTF_VERBOSE_LEVEL=0
export MAX_WORKGROUP_SIZE=2
export DTF_SYNCH_WALLTIME=1
rm $DTF_GLOBAL_PATH/port_*

#mpirun -np 1  valgrind --leak-check=full ./writer &
#mpirun -np 1  valgrind --leak-check=full ./reader

mpirun -np 4 ./writer &
mpirun -np 4 ./reader
