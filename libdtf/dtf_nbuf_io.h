#ifndef DTF_NBUF_IO_H_INCLUDED
#define DTF_NBUF_IO_H_INCLUDED

#include <mpi.h>

MPI_Offset nbuf_read_write_var(file_buffer_t *fbuf,
                               int varid,
                               const MPI_Offset *start,
                               const MPI_Offset *count,
                               const MPI_Offset *stride,
                               const MPI_Offset *imap,
                               MPI_Datatype dtype,
                               void *buf,
                               int rw_flag,
                               int *request);

#endif /*dtf_NBUF_IO_H_INCLUDED*/
