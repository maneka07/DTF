#ifndef FARB_H_INCLUDED
#define FARB_H_INCLUDED

#ifndef _EXTERN_C_
#ifdef __cplusplus
#define _EXTERN_C_ extern "C"
#else /* __cplusplus */
#define _EXTERN_C_
#endif /* __cplusplus */
#endif /* _EXTERN_C_ */

#include <stddef.h>

#ifndef FARB_IO_MODE_UNDEFINED
#define FARB_IO_MODE_UNDEFINED    1
#define FARB_IO_MODE_FILE         2
#define FARB_IO_MODE_MEMORY       3
#endif

#define FARB_READ   1
#define FARB_WRITE  2

/*Interfaces used by an application*/
_EXTERN_C_ int farb_init(const char *filename, char *module_name);
_EXTERN_C_ int farb_finalize();

/*Interfaces to be used by a File I/O library*/
//_EXTERN_C_ size_t farb_write(const char* filename, const off_t offset, const size_t data_sz, void *const data);
//_EXTERN_C_ size_t farb_read(const char* filename, const off_t  offset, const size_t data_sz, void *const data);
_EXTERN_C_ int farb_write_flag(const char* filename);
_EXTERN_C_ int farb_read_flag(const char* filename);
_EXTERN_C_ void farb_progress_io();
_EXTERN_C_ void farb_open(const char* filename);
_EXTERN_C_ void farb_close(const char* filename);
_EXTERN_C_ int farb_io_mode(const char* filename);
_EXTERN_C_ int farb_def_var(const char* filename, int varid, int ndims, MPI_Offset *shape);
_EXTERN_C_ void farb_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header);
_EXTERN_C_ MPI_Offset farb_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
_EXTERN_C_ MPI_Offset farb_read_write_var(const char *filename, int varid, const MPI_Offset *start,const MPI_Offset *count,const MPI_Offset *stride,const MPI_Offset *imap, MPI_Datatype dtype, void *buf, int rw_flag);

/*     Fortran interfaces    */
void farb_init_(const char *filename, char *module_name, int* ierr);
void farb_finalize_(int* ierr);

#endif // FARB_H_INCLUDED
