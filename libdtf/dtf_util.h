#ifndef DTF_UTIL_H_INCLUDED
#define DTF_UTIL_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include <string.h>
#include <assert.h>

#define MAX_WORKGROUP_SIZE     64

#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0

#define MAX_COMP_NAME 32
#define ASCIILINESZ   1024

#define FILE_READY_TAG      0       /*writer rank -> reader rank*/
#define IO_DATA_REQ_TAG     1       /*master -> writer*/
#define IO_DATA_TAG         2       /*writer -> reader*/
#define READ_DONE_TAG       3
#define SKIP_MATCH_TAG   	4      
#define FILE_INFO_TAG       5
#define FILE_INFO_REQ_TAG   6
#define MATCH_DONE_TAG      7
#define IO_REQS_TAG         8
#define COMP_SYNC_TAG	9

#define IODB_BUILD_VARID    0  /*Distribute ioreqs based on var id*/
#define IODB_BUILD_BLOCK    1  /*Distribute ioreqs by dividing var to blocks*/

#define VERBOSE_ERROR_LEVEL   0
#define VERBOSE_DBG_LEVEL     1
#define VERBOSE_ALL_LEVEL     2

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

typedef struct component{
    unsigned int    id;
    char            name[MAX_COMP_NAME];
    int             connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    MPI_Comm        comm;   /*intra or inter component communicator*/
}component_t;

typedef struct dtf_config{
    int         distr_mode;
    int         buffered_req_match;    /*Should we buffer the data if request matching is enabled?*/
    int         data_msg_size_limit;
    int         detect_overlap_flag;    /*should master process detect write overlap by different processes?*/
    int         do_checksum;
    MPI_Offset  iodb_range;  		  /*the size of the data block in the first dimension*/
    int         iodb_build_mode;      /*IODB_BUILD_VARID - based on var ids, IODB_BUILD_RANGE - based on data block range*/
    int         log_ioreqs;
}dtf_config_t;

typedef struct stats{
    int             ndata_msg_sent;
    double          accum_dbuff_time;
    size_t          accum_dbuff_sz;
    double          accum_comm_time;
    double          accum_hdr_time;
    double          parse_ioreq_time;
    size_t          data_msg_sz;      /*Accumulated size of messages */
    double          accum_match_time;  /*Total time spent in match_ioreqs*/
    double          idle_time;
    double          idle_do_match_time;
    double          accum_progr_time;
    double          accum_do_matching_time;
    double          accum_rw_var;
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
    
    double          st_mtch_hist;
    double          end_mtch_hist;
    double          st_mtch_rest;
    double          end_mtch_rest;
    double          st_fin;
} stats_t;

typedef struct dtf_msg{
    MPI_Request req;
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
extern struct file_buffer* gl_filebuf_list;        /*List of all file buffers*/
extern struct fname_pattern *gl_fname_ptrns;    /*Patterns for file name*/
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
extern dtf_msg_t *gl_msg_q;
extern file_info_t *gl_finfo_list;
char _buff[1024];


char error_string[1024];


/*FUNCTIONS*/
int match_ptrn(char* pattern, const char* filename, char** excl_fnames, int nexcls);
void notify_file_ready(file_buffer_t *fbuf);
void close_file(file_buffer_t *fbuf);
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header);
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void open_file(file_buffer_t *fbuf, MPI_Comm comm);
MPI_Offset to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord);
void* dtf_malloc(size_t size);
void dtf_free(void *ptr, size_t size);
void find_fit_block(int ndims,
		    int cur_dim,
		    const MPI_Offset *start,
		    const MPI_Offset *count,
		    MPI_Offset *cur_start,
		    MPI_Offset *cur_count,
		    const size_t sbufsz,
		    const size_t el_sz,
		    MPI_Offset *cur_nelems,
		    MPI_Offset tot_nelems);

void get_put_data(dtf_var_t *var,
                  MPI_Datatype dtype,
                  unsigned char *block_data,
                  const MPI_Offset *block_start,
                  const MPI_Offset *block_count,
                  const MPI_Offset subbl_start[],
                  const MPI_Offset subbl_count[],
                  unsigned char *subbl_data,
                  int get_put_flag,
                  int convert_flag);

void recur_get_put_data(dtf_var_t *var,
                          MPI_Datatype dtype,
                          unsigned char *block_data,
                          const MPI_Offset *block_start,
                          const MPI_Offset *block_count,
                          const MPI_Offset subbl_start[],
                          const MPI_Offset subbl_count[],
                          int dim,
                          MPI_Offset coord[],
                          unsigned char *subbl_data,
                          int get_put_flag,
                          int convert_flag);
void shift_coord(int ndims, const MPI_Offset *bl_start,
                 const MPI_Offset *bl_count, MPI_Offset *subbl_start,
                 MPI_Offset *subbl_count, MPI_Offset fit_nelems);
void convertcpy(MPI_Datatype type1, MPI_Datatype type2, void* srcbuf, void* dstbuf, int nelems);
double compute_checksum(void *arr, int ndims, const MPI_Offset *shape, MPI_Datatype dtype);
dtf_msg_t *new_dtf_msg(void *buf, size_t bufsz, int tag);
void delete_dtf_msg(dtf_msg_t *msg);
void print_stats();
int inquire_root(const char *filename);
#endif
