#ifndef DTF_MEM_H_INCLUDED
#define DTF_MEM_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include "dtf_req_match.h"

//MPI_Offset mem_write(dtf_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_contiguous_write(dtf_var_t *var, MPI_Offset offset, MPI_Offset data_sz, void *data);
MPI_Offset mem_contiguous_read(dtf_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
void get_contig_mem_list(dtf_var_t *var,
                         MPI_Datatype dtype,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         int *nelems,
                         contig_mem_chunk_t **list);
//MPI_Offset mem_read(dtf_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data);
MPI_Offset mem_noncontig_write(dtf_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *buf);
MPI_Offset mem_noncontig_read(dtf_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data);
#endif // dtf_MEM_H_INCLUDED
