#!/bin/bash


export FARB_GLOBAL_PATH=.
export LD_LIBRARY_PATH=../../libfarb:$HOME/apps/pnetcdf/lib:$LD_LIBRARY_PATH
export FARB_VERBOSE_LEVEL=4

rm $FARB_GLOBAL_PATH/port_*

mpirun -np 2 ./writer &
mpirun -np 2 ./reader
