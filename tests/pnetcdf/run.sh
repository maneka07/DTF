#!/bin/bash


export FARB_GLOBAL_PATH=.
export LD_LIBRARY_PATH=../../libfarb:$HOME/apps/pnetcdf/lib:$LD_LIBRARY_PATH
export FARB_VERBOSE_LEVEL=1

rm $FARB_GLOBAL_PATH/port_*
rm restart_*

mpirun -np 2 ./writer &
mpirun -np 2 ./reader
