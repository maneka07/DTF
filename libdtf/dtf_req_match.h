#ifndef DTF_REQ_MATCH_H_INCLUDED
#define DTF_REQ_MATCH_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include "rb_red_black_tree.h"

#define DEFAULT_BLOCK_SZ_RANGE 1
//#define AUTO_BLOCK_SZ_RANGE    0

#define DTF_DATA_MSG_SIZE_LIMIT 256*1024*1024


typedef struct io_req{
    unsigned int            id;
    int                     var_id;
    int                     sent_flag;  /*set to 1 if this req has already been forwarded to the master*/
    int                     rw_flag;
    void                    *user_buf;
    MPI_Datatype            dtype;                /*need to save the type passed in the request
                                                    in case there is mismatch with the type passed
                                                    when the var was defined*/
    MPI_Offset              user_buf_sz;
    MPI_Offset              *start;
    MPI_Offset              *count;
    MPI_Offset              get_sz;         /*size of data received from writer ranks*/
    unsigned                is_buffered;    /*1 if the data is buffered. 0 otherwise*/
    unsigned                is_permanent;  /*needed for SCALE-LETKF as the timeframe var is constantly overwritten in the new iteration but we need to keep the last value*/
    double                  checksum;
    struct io_req           *next;
    struct io_req           *prev;
}io_req_t;

/*All write records will be grouped by var_id.
  All read record will be grouped by reader's rank.
  This is done for simplicity of matching and sending
  data requests to writer ranks*/



typedef struct write_db_item{
    int                    var_id;
    int                    ndims;
    void                   *dblocks; /*rb_red_blk_tree for ndims>0, block_t for ndims=0*/
    MPI_Offset             nblocks;
    struct write_db_item   *next;
}write_db_item_t;


typedef struct read_dblock{
    int                     var_id;
    int                     ndims;
    MPI_Offset              *start;
    MPI_Offset              *count;
    struct read_dblock   *next;
    struct read_dblock   *prev;
}read_dblock_t;

typedef struct read_db_item{
    int                     rank;
    MPI_Comm                comm;   /*Both writer and reader may read the file, hence, need to distinguish.*/
    read_dblock_t           *dblocks;
    read_dblock_t           *last_block;
    MPI_Offset              nblocks;
    struct read_db_item     *next;
    struct read_db_item     *prev;
}read_db_item_t;

typedef struct ioreq_db{
    int                  updated_flag;
    MPI_Offset           nritems;
    //TODO replace list with something more efficient
    struct write_db_item *witems;
    struct read_db_item  *ritems;
}ioreq_db_t;

typedef struct master_info{
    unsigned int         nwranks_completed;  /*Number of writer ranks that completed their read requests
                                              (writer can also read the file). Counted only by master 0.*/
    unsigned int         nwranks_opened;     /*Number of writer ranks that opened the file*/
    struct ioreq_db      *iodb;
    int is_master_flag;  /*is this rank a master rank*/
    int nmasters;   /*Number of master nodes that hold data for request matching*/
    int *masters;   /*Ranks of master nodes on the writer's side*/
    unsigned int my_wg_sz;
    int *my_wg;
} master_info_t;

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
                    int var_id,
                    int ndims,
                    MPI_Datatype dtyp,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    void *buf,
                    int rw_flag,
                    int buffered);



void progress_msg_queue();
void add_ioreq(io_req_t **list, io_req_t *ioreq);
//void delete_ioreqs(file_buffer_t *fbuf);
void delete_ioreqs(file_buffer_t *fbuf, int finalize);
void progress_io_matching();
void send_ioreqs(file_buffer_t *fbuf, int intracomp_match);
void clean_iodb(ioreq_db_t *iodb);
int  match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag);
//void match_ioreqs_all(int rw_flag);
int  init_req_match_masters(MPI_Comm comm, master_info_t *mst_info);
void init_iodb(file_buffer_t *fbuf);
void unpack_file_info(MPI_Offset bufsz, void *buf);
void send_file_info(file_buffer_t *fbuf, int reader_root);
void notify_complete_multiple(file_buffer_t *fbuf);
#endif // dtf_REQ_MATCH_H_INCLUDED
