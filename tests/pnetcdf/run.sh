#!/bin/bash

export FARB_GLOBAL_PATH=.
export LD_LIBRARY_PATH=../../libfarb:$HOME/apps/pnetcdf/lib:$LD_LIBRARY_PATH
export FARB_VERBOSE_LEVEL=0
export MAX_WORKGROUP_SIZE=2
rm $FARB_GLOBAL_PATH/port_*

#mpirun -np 1  valgrind --leak-check=full ./writer &
#mpirun -np 1  valgrind --leak-check=full ./reader

mpirun -np 16  ./writer &
mpirun -np 16 ./reader
