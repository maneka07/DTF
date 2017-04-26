#ifndef DTF_UTIL_H_INCLUDED
#define DTF_UTIL_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"

#define MAX_WORKGROUP_SIZE     64

#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0

#define MAX_COMP_NAME 32
#define ASCIILINESZ   1024

#define FILE_READY_TAG      0       /*writer rank -> reader rank*/
#define IO_DATA_REQ_TAG     1       /*master -> writer*/
#define DONE_MULTIPLE_FLAG  2
#define IO_DATA_TAG         3       /*writer -> reader*/
#define READ_DONE_TAG       4       /*reader->master*/
#define IO_CLOSE_FILE_TAG   5      /*reader->master, master->writers*/
#define FILE_INFO_TAG       6
#define FILE_INFO_REQ_TAG   7
#define ROOT_MST_TAG        8
#define MATCH_DONE_TAG      9
#define IO_REQS_TAG         10

#define DISTR_MODE_REQ_MATCH    1

typedef struct component{
    unsigned int    id;
    char            name[MAX_COMP_NAME];
    int             connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    MPI_Comm        comm;   /*intra or inter component communicator*/
}component_t;

typedef struct dtf_config{
    int         distr_mode;
    int         buffered_req_match;    /*Should we buffer the data if request matching is enabled?*/
    int         io_db_type;  /*data blocks(0)*/
    int         data_msg_size_limit;
    int         detect_overlap_flag;    /*should master process detect write overlap by different processes?*/
}dtf_config_t;

typedef struct stats{
    int             nmsg_sent;
    int             nmatching_msg_sent;
    double          accum_comm_time;
    size_t          accum_msg_sz;      /*Accumulated size of messages */
    double          accum_match_time;  /*Total time spent in match_ioreqs*/
    double          accum_db_match_time; /*Time that master processes spend matching reqs from iodb*/
    int             ndb_match;       /*How many times I/O req matching was performed*/
    unsigned int    num_tsrch;
    double          t_treesrch;
    double          t_progress;
    double          walltime;
    size_t          malloc_size;
    unsigned long   nprogress_call;
    unsigned        nioreqs;
} stats_t;


struct file_buffer;
int mpitype2int(MPI_Datatype dtype);
MPI_Datatype int2mpitype(int num);

void notify_file_ready(file_buffer_t *fbuf);
void close_file(file_buffer_t *fbuf);
int file_buffer_ready(const char* filename);
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header);
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
int set_distr_count(file_buffer_t *fbuf, int varid, int count[]);
void open_file(file_buffer_t *fbuf, MPI_Comm comm);
MPI_Offset to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord);
void* dtf_malloc(size_t size);
void dtf_free(void *ptr, size_t size);
void process_file_info_req_queue();
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

#endif
