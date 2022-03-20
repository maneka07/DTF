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

#ifndef DTF_VAR_H_INCLUDED
#define DTF_VAR_H_INCLUDED

#include <mpi.h>

struct file_buffer;
struct io_req;

typedef struct dtf_var{
    int                     id;         /* varid assigned by pnetcdf*/
    MPI_Datatype            dtype;      /*Datatype of the variable*/
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
    struct io_req           *ioreqs;    /*Read or write I/O requests*/
    double                  checksum;   /*to check if what was written to the var is what was read*/
 }dtf_var_t;

dtf_var_t* 	new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void 		add_var(struct file_buffer *fbuf, dtf_var_t *var);
MPI_Offset 	read_write_var(struct file_buffer *fbuf,
								   int varid,
								   const MPI_Offset *start,
								   const MPI_Offset *count,
								   MPI_Datatype dtype,
								   void *buf,
								   int rw_flag);
int 		def_var(struct file_buffer *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void 		delete_var(struct file_buffer *fbuf, dtf_var_t* var);

#endif
