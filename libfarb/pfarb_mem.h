#ifndef PFARB_MEM_H_INCLUDED
#define PFARB_MEM_H_INCLUDED

#include <mpi.h>
#include "pfarb_file_buffer.h"

//MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_contiguous_write(farb_var_t *var, MPI_Offset offset, MPI_Offset data_sz, void *data);
MPI_Offset mem_contiguous_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
void get_contig_mem_list(farb_var_t *var,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         int *nelems,
                         contig_mem_chunk_t **list);
//MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_noncontig_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *buf);
MPI_Offset mem_noncontig_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data);
#endif // PFARB_MEM_H_INCLUDED
