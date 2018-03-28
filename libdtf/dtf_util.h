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
#define READ_DONE_CONFIRM_TAG       4
#define FILE_INFO_TAG       5
#define FILE_INFO_REQ_TAG   6
#define MATCH_DONE_TAG      7
#define IO_REQS_TAG         8
#define COMP_SYNC_TAG		9
#define COMP_FINALIZED_TAG   10

#define IODB_BUILD_VARID    0  /*Distribute ioreqs based on var id*/
#define IODB_BUILD_BLOCK    1  /*Distribute ioreqs by dividing var to blocks*/

#define VERBOSE_ERROR_LEVEL   0
#define VERBOSE_DBG_LEVEL     1
#define VERBOSE_ALL_LEVEL     2
#define VERBOSE_RB_TREE_LEVEL 3

#define DTF_TIMEOUT       60 

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
    if(gl_verbose >= dbg_level){  \
                memset(_buff,0,1024);                         \
                snprintf(_buff,1024,__VA_ARGS__);             \
                fprintf(stdout, "%s %d [%.3f]: %s\n", gl_my_comp_name, gl_my_rank, MPI_Wtime() - gl_stats.walltime, _buff);  \
                fflush(stdout);   \
    }           \
}while(0)

#define CHECK_MPI(errcode) do{   \
        if (errcode != MPI_SUCCESS) {   \
           int length_of_error_string;  \
           MPI_Error_string(errcode, error_string, &length_of_error_string);  \
           DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error is: %s", error_string);       \
           MPI_Abort(MPI_COMM_WORLD, errcode);                              \
        }                                                                   \
} while(0)



typedef struct dtf_config{
    int         distr_mode;
    int         buffer_data;    /*Should we buffer the data if request matching is enabled?*/
    int         data_msg_size_limit;
    int         detect_overlap_flag;    /*should master process detect write overlap by different processes?*/
    int         do_checksum;
    MPI_Offset  iodb_range;  		  /*the size of the data block in the first dimension*/
    int         iodb_build_mode;      /*IODB_BUILD_VARID - based on var ids, IODB_BUILD_RANGE - based on data block range*/
    int         log_ioreqs;
	double          t_send_ioreqs_freq;  /*time out upon reaching which I/O reqs must be sent*/
}dtf_config_t;

typedef struct stats{
    int             ndata_msg_sent;
    double          accum_dbuff_time;
    size_t          accum_dbuff_sz;
    double          t_comm;
    double          t_hdr;
    double          parse_ioreq_time;
    size_t          data_msg_sz;      /*Accumulated size of messages */
    double          idle_time;
    double          idle_do_match_time;
    double          t_progr_comm;
    double          t_do_match;
    double          t_rw_var;
    double          master_time;  /*measure accum time for master-related work*/
    int             ndb_match;       /*How many times I/O req matching was performed*/
    double          walltime;
    size_t          malloc_size;
    unsigned long   nprogress_call;
    unsigned long   nioreqs;
    unsigned long   iodb_nioreqs;   /*number of blocks in iodb for matching*/
    unsigned        nbl;    /*number of blocks transfered*/
    unsigned        ngetputcall;  /*how many times had to use a subblock extraction function*/
    double          timer_start;   
    double          timer_accum;
    double          user_timer_start;
    double          user_timer_accum;
    int             nfiles;
    
    double          transfer_time;  /*data transfer time=I/O calls+dtf_transfer */
    double          dtf_time;       /*Total time spent inside DTF*/
    double          st_mtch_hist;
    double          end_mtch_hist;
    double          st_mtch_rest;
    double          end_mtch_rest;
    double          st_fin;
    double          t_open_hist;
    double          t_open_rest;
    double          t_mtch_hist;
    double          t_mtch_rest;
} stats_t;

typedef struct dtf_msg{
    MPI_Request req;
    int src;   //if src is DTF_UNDEFINED, this is an outgoing message, otherwise incoming
    void *buf;
    int tag;
    size_t bufsz;
    struct dtf_msg *next;
    struct dtf_msg *prev;
}dtf_msg_t;

typedef struct file_info{
	char filename[MAX_FILE_NAME];
	int root_writer;
	struct file_info *next;
	struct file_info *prev;
}file_info_t;

struct file_buffer;
struct file_info_req_q;
int mpitype2int(MPI_Datatype dtype);
MPI_Datatype int2mpitype(int num);

/*GLOBAL VARIABLES*/
extern file_buffer_t* 	gl_filebuf_list;        /*List of all file buffers*/
extern struct fname_pattern*	gl_fname_ptrns;    /*Patterns for file name*/
extern struct component *gl_comps;                 /*List of components*/
extern int gl_my_comp_id;                          /*Id of this compoinent*/
extern int gl_ncomp;                               /*Number of components*/
extern int gl_verbose;
extern int gl_my_rank;                         /*For debug messages*/
extern int gl_scale;                           /*special flag for scale-letkf execution*/
extern struct dtf_config gl_conf;                 /*Framework settings*/
extern struct stats gl_stats;
extern char *gl_my_comp_name;
extern void* gl_msg_buf;
extern struct file_info_req_q *gl_finfo_req_q;    
extern file_info_t *gl_finfo_list;
char _buff[1024];

char error_string[1024];


/*FUNCTIONS*/
void 		notify_file_ready(file_buffer_t *fbuf);
void 		send_mst_info(file_buffer_t *fbuf, int tgt_root, int tgt_comp);
void* 		dtf_malloc(size_t size);
void  		dtf_free(void *ptr, size_t size);
void  		convertcpy(MPI_Datatype type1, MPI_Datatype type2, void* srcbuf, void* dstbuf, int nelems);
double 		compute_checksum(void *arr, int ndims, const MPI_Offset *shape, MPI_Datatype dtype);
dtf_msg_t*	new_dtf_msg(void *buf, size_t bufsz, int src, int tag);
void 		print_stats();
int 		inquire_root(const char *filename);
MPI_Offset 	to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord);
int 		boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start,const MPI_Offset *count );
void 		get_put_data(dtf_var_t *var,
                  MPI_Datatype dtype,
                  unsigned char *block_data,
                  const MPI_Offset *block_start,
                  const MPI_Offset *block_count,
                  const MPI_Offset subbl_start[],
                  const MPI_Offset subbl_count[],
                  unsigned char *subbl_data,
                  int get_put_flag,
                  int convert_flag);
void progress_send_queue();
void progress_recv_queue();

#endif
