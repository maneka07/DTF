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

#ifndef DTF_REQ_MATCH_H_INCLUDED
#define DTF_REQ_MATCH_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include "rb_red_black_tree.h"

#define DEFAULT_BLOCK_SZ_RANGE 1
//#define AUTO_BLOCK_SZ_RANGE    0

#define DTF_DATA_MSG_SIZE_LIMIT 64*1024*1024

#define DTF_GROUP_MST 0
#define DTF_GROUP_WG  1


/*This structure is used for debugging purposes:
 * log ioreqs when file i/o is used,
 * then compute the checksum of the user buffer*/
typedef struct io_req_log{
	int				var_id;
	int 			ndims;
	MPI_Datatype 	dtype;
	int				rw_flag;
	void 			*user_buf;
	MPI_Offset      req_data_sz;
	MPI_Offset 		*start;
	MPI_Offset 		*count;
	MPI_Offset      *glob_arr_size;  /*In case dtype is a subarray (derived datatype), we need to know the size of the global array to copy the data correctly*/
	struct io_req_log *next;
}io_req_log_t;

typedef struct dtype_params{
	MPI_Offset  *orig_array_size;
	MPI_Offset  *orig_start; 
} dtype_params_t;

typedef struct io_req{
    unsigned long           id;  
    int                     sent_flag;  /*set to 1 if this req has already been forwarded to the master*/
    int                     rw_flag;
    void                    *user_buf;
    MPI_Datatype            dtype;                /*save the type passed in the request
                                                    in case there is mismatch with the type passed
                                                    when the var was defined*/
    MPI_Offset              req_data_sz;
    MPI_Offset              *start;
    MPI_Offset              *count;
    dtype_params_t          *derived_params;  /*In case it's a derived datatype, need to know original array size and shift within the array where the data should be written to/read from*/
    MPI_Offset              get_sz;         /*size of data received from writer ranks*/
    unsigned                is_buffered;    /*1 if the data is buffered. 0 otherwise*/
    double                  checksum;
    struct io_req           *next;
    struct io_req           *prev;
}io_req_t;

/*All write records will be grouped by var_id.
  All read record will be grouped by reader's rank.
  This is done for simplicity of matching and sending
  data requests to writer ranks*/

typedef struct file_info_req_q{
    char filename[MAX_FILE_NAME];
    void *buf;  /*consists of
                  - filename [char[]]
                  - root reader rank to whom file info should be sent [int]
                 */
    struct file_info_req_q *next;
    struct file_info_req_q *prev;
}file_info_req_q_t;

io_req_t *new_ioreq(int id,
                    int ndims,
                    MPI_Datatype dtype,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    dtype_params_t *derived_params,
                    void *buf,
                    int rw_flag,
                    int buffered);
void add_ioreq(io_req_t **list, io_req_t *ioreq);
void delete_ioreqs(file_buffer_t *fbuf);
void progress_comm(int ignore_idle);
int  match_ioreqs(file_buffer_t *fbuf);
void match_ioreqs_all_files();
void send_file_info(file_buffer_t *fbuf, int reader_root);
void notify_complete_multiple(file_buffer_t *fbuf);
void log_ioreq(file_buffer_t *fbuf,
			  int varid, int ndims,
			  const MPI_Offset *start,
			  const MPI_Offset *count,
			  MPI_Datatype dtype,
			  void *buf,
			  int rw_flag);
void send_data_blocking(file_buffer_t *fbuf, void* buf, int bufsz);
void send_data_nonblocking(file_buffer_t *fbuf, void* buf, int bufsz);
int  parse_msg(int comp, int src, int tag, void *rbuf, int bufsz, int is_queued);
void send_ioreqs_by_block(file_buffer_t *fbuf);
void send_ioreqs_by_var(file_buffer_t *fbuf);
#endif // dtf_REQ_MATCH_H_INCLUDED
