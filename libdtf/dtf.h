#ifndef DTF_H_INCLUDED
#define DTF_H_INCLUDED

#ifndef _EXTERN_C_
#ifdef __cplusplus
#define _EXTERN_C_ extern "C"
#else /* __cplusplus */
#define _EXTERN_C_
#endif /* __cplusplus */
#endif /* _EXTERN_C_ */

#define DTF_IO_MODE_FILE         1
#define DTF_IO_MODE_MEMORY       2

#define DTF_READ   1
#define DTF_WRITE  2

#define DTF_UNDEFINED      -1

#define DTF_UNLIMITED  0L  /*Unlimited dimension*/

/*User API*/
_EXTERN_C_ int dtf_init(const char *filename, char *module_name);
_EXTERN_C_ int dtf_finalize();
_EXTERN_C_ void dtf_print(const char *str);
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid);
_EXTERN_C_ void dtf_match_multiple(int ncid);
_EXTERN_C_ void dtf_tstart();
_EXTERN_C_ void dtf_tend();
_EXTERN_C_ void dtf_time_start();
_EXTERN_C_ void dtf_time_tend();
_EXTERN_C_ int  dtf_transfer(const char *filename, int ncid );
_EXTERN_C_ int  dtf_transfer_v2(const char *filename, int ncid, int it );

/*     Fortran interfaces    */
void dtf_tend_();
void dtf_tstart_();
void dtf_time_end_();
void dtf_time_start_();
void dtf_init_(const char *filename, char *module_name, int* ierr);
void dtf_finalize_(int* ierr);
void dtf_transfer_(const char *filename, int *ncid, int *ierr);
void dtf_print_(const char *str);
void dtf_print_data_(int *varid, int *dtype, int *ndims, MPI_Offset* count, void* data);

/*Interfaces used by PnetCDF*/
_EXTERN_C_ void dtf_open(const char* filename, int omode, MPI_Comm comm);
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
                                          int rw_flag);
_EXTERN_C_ void dtf_log_ioreq(const char *filename,
                                          int varid, int ndims,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag);
_EXTERN_C_ void dtf_enddef(const char *filename);
_EXTERN_C_ void dtf_print_data(int varid, int dtype, int ndims, MPI_Offset* count, void* data);
_EXTERN_C_ void dtf_set_ncid(const char *filename, int ncid);

#endif // DTF_H_INCLUDED
