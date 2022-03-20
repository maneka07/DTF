/*
 © Copyright 2019 RIKEN Center for Computational Science, System Software
 *	Development Team All rights reserved.

 * The Data Transfer Framework (DTF) was designed to work with 
 * the PnetCDF library developed by Northwestern University and 
 * Argonne National Laboratory. The modified version of PnetCDF is 
 * provided with this distribution.
 * 
 * Access and use of this software shall impose the following obligations 
 * and understandings on the user. The user is granted the right, without 
 * any fee or cost, to use, copy, modify, alter, enhance and distribute 
 * this software, and any derivative works thereof, and its supporting 
 * documentation for any purpose whatsoever, provided that this entire 
 * notice appears in all copies of the software, derivative works and 
 * supporting documentation.  Further, RIKEN requests that the user credit 
 * RIKEN in any publications that result from the use of this software or 
 * in any product that includes this software.  The name RIKEN, however, 
 * may not be used in any advertising or publicity to endorse or promote 
 * any products or commercial entity unless specific written permission is 
 * obtained from RIKEN. The user also understands that RIKEN is not 
 * obligated to provide the user with any support, consulting, training or 
 * assistance of any kind with regard to the use, operation and 
 * performance of this software nor to provide the user with any updates, 
 * revisions, new versions or "bug fixes."
 * 
 * THIS SOFTWARE IS PROVIDED BY RIKEN "AS IS" AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN 
 * NO EVENT SHALL RIKEN BE LIABLE FOR ANY SPECIAL, INDIRECT OR 
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF 
 * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR 
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, 
 * USE OR PERFORMANCE OF THIS SOFTWARE.
 *
*/

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

/*		User API		*/
_EXTERN_C_ int dtf_init(const char *filename, char *module_name);
_EXTERN_C_ int dtf_finalize();
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid);
_EXTERN_C_ void dtf_transfer_multiple(const char *filename, int ncid);
_EXTERN_C_ void dtf_tstart();
_EXTERN_C_ void dtf_tend();
_EXTERN_C_ void dtf_time_start();
_EXTERN_C_ void dtf_time_tend();
_EXTERN_C_ int  dtf_transfer(const char *filename, int ncid );
_EXTERN_C_ int  dtf_transfer_all_files();


/*     Fortran interfaces    */
void dtf_time_end_();
void dtf_time_start_();
void dtf_init_(const char *filename, char *module_name, int* ierr);
void dtf_finalize_(int* ierr);
void dtf_transfer_(const char *filename, int *ncid, int *ierr);
void dtf_transfer_all_files_();
void dtf_transfer_multiple_(const char *filename, int *ncid);

/*		Interfaces used by PnetCDF		*/
_EXTERN_C_ void dtf_open(const char* filename, int omode, MPI_Comm comm);
_EXTERN_C_ void dtf_close(const char* filename);
_EXTERN_C_ void dtf_create(const char *filename, MPI_Comm comm);
_EXTERN_C_ int  dtf_io_mode(const char* filename);
_EXTERN_C_ int  dtf_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
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
_EXTERN_C_ void dtf_set_ncid(const char *filename, int ncid);
_EXTERN_C_ MPI_File *dtf_get_tmpfile();

#endif // DTF_H_INCLUDED
