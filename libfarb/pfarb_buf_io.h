#ifndef PFARB_BUF_IO_H_INCLUDED
#define PFARB_BUF_IO_H_INCLUDED

#include "pfarb_file_buffer.h"
#include <mpi.h>

//void pack_vars(file_buffer_t *fbuf, int dst_rank, int *buf_sz, void **buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord);
//void unpack_vars(file_buffer_t *fbuf, int buf_sz, void *buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord);
int receive_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm);
int send_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm);

MPI_Offset buf_read_write_var( file_buffer_t *fbuf,
                               int varid,
                               const MPI_Offset *start,
                               const MPI_Offset *count,
                               const MPI_Offset *stride,
                               const MPI_Offset *imap,
                               MPI_Datatype dtype,
                               void *buf,
                               int rw_flag);

#endif /*PFARB_BUF_IO_H_INCLUDED*/
