#ifndef DTF_H_INCLUDED
#define DTF_H_INCLUDED

#ifndef _EXTERN_C_
#ifdef __cplusplus
#define _EXTERN_C_ extern "C"
#else /* __cplusplus */
#define _EXTERN_C_
#endif /* __cplusplus */
#endif /* _EXTERN_C_ */


#ifndef DTF_IO_MODE_UNDEFINED
#define DTF_IO_MODE_UNDEFINED    1
#define DTF_IO_MODE_FILE         2
#define DTF_IO_MODE_MEMORY       3
#endif

#define DTF_READ   1
#define DTF_WRITE  2

#define DTF_UNLIMITED  0L  /*Unlimited dimension*/

/*TODO: implement hierarchical cleaning in case error occurs*/

/*Interfaces used by an application*/
_EXTERN_C_ int dtf_init(const char *filename, char *module_name);
_EXTERN_C_ int dtf_finalize();
_EXTERN_C_ int dtf_set_distr_count(const char* filename, int varid, int count[]);
_EXTERN_C_ int dtf_wait_all(const char* filename);
_EXTERN_C_ int dtf_wait(const char* filename, int request);
/*Interfaces to be used by a File I/O library*/
_EXTERN_C_ int dtf_write_flag(const char* filename);
_EXTERN_C_ int dtf_read_flag(const char* filename);
_EXTERN_C_ void dtf_open(const char* filename, MPI_Comm comm);
_EXTERN_C_ void dtf_close(const char* filename);
_EXTERN_C_ void dtf_create(const char *filename, MPI_Comm comm, int ncid);
_EXTERN_C_ int dtf_io_mode(const char* filename);
_EXTERN_C_ int dtf_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
_EXTERN_C_ void dtf_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header);
_EXTERN_C_ MPI_Offset dtf_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
_EXTERN_C_ MPI_Offset dtf_read_write_var(const char *filename,
                                          int varid,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          const MPI_Offset *stride,
                                          const MPI_Offset *imap,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag,
                                          int *request);
_EXTERN_C_ int dtf_match_ioreqs(const char* filename);
_EXTERN_C_ int dtf_match_io(const char *filename, int ncid, int intracomp_io_flag);
_EXTERN_C_ void dtf_match_io_all(int rw_flag);
_EXTERN_C_ void dtf_enddef(const char *filename);

/*     Fortran interfaces    */
void dtf_init_(const char *filename, char *module_name, int* ierr);
void dtf_finalize_(int* ierr);
void dtf_match_io_(const char *filename, int *ncid, int *intracomp_io_flag, int *ierr);
void dtf_match_io_all_(int *rw_flag);


#endif // DTF_H_INCLUDED
