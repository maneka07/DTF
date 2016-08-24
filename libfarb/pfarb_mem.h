#ifndef PFARB_MEM_H_INCLUDED
#define PFARB_MEM_H_INCLUDED

#include "pfarb_util.h"

MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);


#endif // PFARB_MEM_H_INCLUDED
