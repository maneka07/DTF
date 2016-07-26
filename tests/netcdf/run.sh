#!/bin/bash

NETCDF_PATH=
export FARB_GLOBAL_PATH=
export LD_LIBRARY_PATH=../../comm_lib/libfarb:$HOME/apps/netcdf/lib:$LD_LIBRARY_PATH
export FARB_VERBOSE_LEVEL=1

rm $FARB_GLOBAL_PATH/port_*
mpirun -np 2 ./writer &
mpirun -np 2 ./reader
