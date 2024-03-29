﻿/*
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

#ifndef DTF_UTIL_H_INCLUDED
#define DTF_UTIL_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include "dtf_component.h"
#include <string.h>
#include <assert.h>

#define MAX_WORKGROUP_SIZE     64

#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0
#define ASCIILINESZ   1024

#define FILE_READY_TAG      0       /*writer rank -> reader rank*/
#define IO_DATA_REQ_TAG     1       /*master -> writer*/
#define IO_DATA_TAG         2       /*writer -> reader*/
#define READ_DONE_TAG       3
#define DONE_MULTIPLE_TAG   4
#define FILE_INFO_TAG       5
#define FILE_INFO_REQ_TAG   6
#define MATCH_DONE_TAG      7
#define IO_REQS_TAG         8
#define COMP_FINALIZED_TAG  9
#define COMP_ABORT_TAG      10

#define IODB_BUILD_VARID    0  /*Distribute ioreqs based on var id*/
#define IODB_BUILD_BLOCK    1  /*Distribute ioreqs by dividing var to blocks*/

#define VERBOSE_ERROR_LEVEL   0
#define VERBOSE_DBG_LEVEL     1
#define VERBOSE_ALL_LEVEL     2
#define VERBOSE_RB_TREE_LEVEL 3
#define VERBOSE_OVERRIDE_LEVEL 4
/*NOTE: These two definitions are copied from pnetcdf.h
 * since I wanted to be able to compile DTF without having 
 * to link it to pnetcdf.*/
#define NC_NOWRITE	 0x0000	/**< Set read-only access for nc_open(). */
#define NC_WRITE    	 0x0001	/**< Set read-write access for nc_open(). */

#define ENQUEUE_ITEM(item, queue) do{\
    if(queue == NULL)   \
        queue = item;    \
    else{   \
        item->next = queue;  \
        queue->prev = item;  \
        queue = item;    \
    }   \
    DTF_DBG(VERBOSE_DBG_LEVEL, "enq_item %p", (void*)item);    \
} while(0)

#define DEQUEUE_ITEM(item, queue) do{   \
    if(item->prev != NULL)  \
        item->prev->next = item->next;  \
    if(item->next != NULL)  \
        item->next->prev = item->prev;  \
    if(queue == item){   \
        queue = item->next; \
    }   \
    DTF_DBG(VERBOSE_DBG_LEVEL, "deq_item %p", (void*)item);    \
    if(queue == NULL)   \
            DTF_DBG(VERBOSE_DBG_LEVEL, "queue empty");  \
} while(0)



#define DTF_DBG(dbg_level, ...) do{  \
    if(gl_verbose >= dbg_level || dbg_level==VERBOSE_OVERRIDE_LEVEL){  \
                memset(_buff,0,1024);                         \
                snprintf(_buff,1024,__VA_ARGS__);             \
                fprintf(stdout, "%s %d [%.3f]: %s\n", _comp_name, gl_proc.myrank, MPI_Wtime() - gl_proc.walltime, _buff);  \
                fflush(stdout);   \
    }           \
}while(0)

#define CHECK_MPI(errcode) do{   \
        if (errcode != MPI_SUCCESS) {   \
           int length_of_error_string;  \
           MPI_Error_string(errcode, _error_string, &length_of_error_string);  \
           DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error is: %s", _error_string);       \
           MPI_Abort(MPI_COMM_WORLD, errcode);                              \
        }                                                                   \
} while(0)



typedef struct dtf_config{
    int         single_mpirun_mode;
    int         buffer_data;    /*Should we buffer the data if request matching is enabled?*/
    int         data_msg_size_limit;
    int         use_msg_buffer;
    int         do_checksum;
    int 		timeout;        /*How long to wait before consider that DTF stalled?*/
    int         ignore_idle;    /*don't abort process if has ben idle for timeout seconds*/
    MPI_Offset  iodb_range;  		  /*the size of the data block in the first dimension*/
    int         iodb_build_mode;      /*IODB_BUILD_VARID - based on var ids, IODB_BUILD_RANGE - based on data block range*/
    int         log_ioreqs;
}dtf_config_t;

typedef struct stats{
	int             num_fpats;
    int             ndata_msg_sent;
    double          accum_dbuff_time;
    size_t          accum_dbuff_sz;
    double          t_comm;
    size_t          data_msg_sz;      /*Accumulated size of messages */
    double          idle_time;
    double          idle_do_match_time;
    double          t_progr_comm;
    double          t_do_match;
    double			t_parse;
    double			t_search;
    unsigned		nsearch;	
    unsigned        ndomatch;	
    double          master_time;  /*measure accum time for master-related work*/ 
    size_t          malloc_size;
    unsigned long   nioreqs;  /*Number of I/O requests issued by the process*/
    unsigned long   nrreqs;   /*number of read I/O requests in iodb for matching*/
    unsigned long   nwreqs;   /*number of write I/O requests in iodb for matching*/
    double          timer_start;   
    double          timer_accum;
    double          user_timer_start;
    double          user_timer_accum;
    unsigned int    nfiles;
    
    unsigned int    nmasters; //remember number of masters of the last opened file
    
    double          transfer_time;  /*data transfer time=I/O calls+dtf_transfer */
    double          dtf_time;       /*Total time spent inside DTF*/
    double          t_idle;        /*used to figure out when program aborted*/
} stats_t;

typedef struct dtf_msg{
    MPI_Request*	reqs;
    int         	nreqs; 
    int 			src;   //if src is DTF_UNDEFINED, this is an outgoing message, otherwise incoming
    void*			buf;
    int 			tag;
    size_t 			bufsz;
    struct dtf_msg* next;
    struct dtf_msg* prev;
}dtf_msg_t;

typedef struct dtf_proc{
	file_buffer_t* 		filebuf_list;   /*List of all file buffers*/
	fname_pattern_t*	fname_ptrns;    /*Patterns for file name*/
	dtf_config_t		conf;        /*Framework settings*/
    stats_t 			stats_info;       /*Execution statistics*/
    component_t*		comps;                 /*List of components*/
	int 				my_comp;                          /*Id of this compoinent*/
	int 				ncomps;                               /*Number of components*/
	int 				myrank;                         /*For debug messages*/
	char 				tmp_file_path[L_tmpnam];
	MPI_File 			tmpfile;
	void* 				msgbuf;
	double          	walltime; 
}dtf_proc_t;

/*GLOBAL VARIABLES*/
extern int gl_verbose;
extern struct dtf_proc gl_proc;

char _error_string[1024];
char _buff[1024];
char _comp_name[MAX_COMP_NAME];

/*FUNCTIONS*/
void 		notify_file_ready(file_buffer_t *fbuf);
void 		send_mst_info(file_buffer_t *fbuf, int tgt_root, int tgt_comp);
void* 		dtf_malloc(size_t size);
void  		dtf_free(void *ptr, size_t size);
void  		convertcpy(MPI_Datatype type1, MPI_Datatype type2, void* srcbuf, void* dstbuf, int nelems);
double 		compute_checksum(void *arr, int ndims, const MPI_Offset *shape, MPI_Datatype dtype);
dtf_msg_t*	new_dtf_msg(void *buf, size_t bufsz, int src, int tag, int nreqs);
void 		delete_dtf_msg(dtf_msg_t *msg);
void 		print_stats();
int 		inquire_root(const char *filename);
MPI_Offset 	to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord);
int 		boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start,const MPI_Offset *count );
void get_put_data(int ndims, 
				  MPI_Datatype orig_dtype,
                  MPI_Datatype tgt_dtype,
				  const MPI_Offset  block_count[], 
                  unsigned char    *block_data,
                  const MPI_Offset subbl_start[],
                  const MPI_Offset subbl_count[],
                  unsigned char *subbl_data,
                  int get_put_flag,
                  int convert_flag);
void progress_send_queue();
void progress_recv_queue();
void translate_ranks(int *from_ranks,  int nranks, MPI_Comm from_comm, MPI_Comm to_comm, int *to_ranks);
int mpitype2int(MPI_Datatype dtype);
MPI_Datatype int2mpitype(int num);
const char *mpi_type_name(MPI_Datatype dtype);
void create_tmp_file();
void delete_tmp_file();
#endif
