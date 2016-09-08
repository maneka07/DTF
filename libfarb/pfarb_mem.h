#ifndef PFARB_MEM_H_INCLUDED
#define PFARB_MEM_H_INCLUDED

#include "pfarb_util.h"

MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_recursive_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *buf, int dim, MPI_Offset coord[]);
MPI_Offset mem_recursive_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data, int dim, MPI_Offset coord[]);
#endif // PFARB_MEM_H_INCLUDED
