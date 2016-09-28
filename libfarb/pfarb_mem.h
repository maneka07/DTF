#ifndef PFARB_MEM_H_INCLUDED
#define PFARB_MEM_H_INCLUDED

#include "pfarb_util.h"

//MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_write(farb_var_t *var, MPI_Offset first_el_coord[],  MPI_Offset data_sz, void *data);
MPI_Offset mem_read(farb_var_t *var, MPI_Offset first_el_coord[],  MPI_Offset data_sz, void *data);
//MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_noncontig_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *buf);
MPI_Offset mem_noncontig_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data);
#endif // PFARB_MEM_H_INCLUDED
