#ifndef DTF_REQ_MATCH_H_INCLUDED
#define DTF_REQ_MATCH_H_INCLUDED

#include <mpi.h>
#include "dtf_file_buffer.h"
#include "rb_red_black_tree.h"

#define UNLIM_NELEMS_RANGE 256

#define DTF_DB_BLOCKS  0  /*Match based on multidimensional data blocks*/
#define DTF_DB_CHUNKS  1  /*Match based on contigious memory chunks*/

typedef struct contig_mem_chunk{
    MPI_Offset              offset;
    MPI_Offset              usrbuf_offset;
    MPI_Offset              data_sz;
    struct contig_mem_chunk *next;
}contig_mem_chunk_t;

typedef struct io_req{
    unsigned int            id;
    int                     var_id;
    int                     sent_flag;  /*set to 1 if this req has already been forwarded to the master*/
    int                     rw_flag;
    void                    *user_buf;
    //MPI_Offset              user_buf_el_sz;     //In case this does
    MPI_Datatype            dtype;                /*need to save the type passed in the request
                                                    in case there is mismatch with the type passed
                                                    when the var was defined*/
    MPI_Offset              user_buf_sz;
    MPI_Offset              *start;
    MPI_Offset              *count;
    MPI_Offset              get_sz;     /*size of data received from writer ranks*/
    int                     nchunks;
    struct contig_mem_chunk *mem_chunks;
    struct io_req           *next;
    struct io_req           *prev;
}io_req_t;

/*All write records will be grouped by var_id.
  All read record will be grouped by reader's rank.
  This is done for simplicity of matching and sending
  data requests to writer ranks*/

//typedef struct chunk_tuple{
//    MPI_Offset             offset;
//    MPI_Offset             data_sz;
//}chunk_tuple_t;

typedef struct write_chunk_rec{
    int                    rank;
    MPI_Offset             offset;
    MPI_Offset             data_sz;
}write_chunk_rec_t;

typedef struct write_dblock{
    int rank;
    MPI_Offset *start;
    MPI_Offset *count;
    struct write_dblock *next;
}write_dblock_t;

typedef struct write_db_item{
    int                    var_id;
    int                    ndims;
    rb_red_blk_tree        *chunks;
    MPI_Offset             nchunks;
    write_dblock_t         *dblocks;
    write_dblock_t         *last_block;
    MPI_Offset             nblocks;
    //struct write_chunk_rec *chunks;
    //struct write_chunk_rec *last;
    struct write_db_item   *next;
}write_db_item_t;

typedef struct read_chunk_rec{
    int                     var_id;
    MPI_Offset             offset;
    MPI_Offset             data_sz;
    struct read_chunk_rec   *next;
    struct read_chunk_rec   *prev;
}read_chunk_rec_t;

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
    MPI_Offset              nchunks;
    struct read_chunk_rec   *chunks;
    struct read_chunk_rec   *last;

    read_dblock_t           *dblocks;
    read_dblock_t           *last_block;
    MPI_Offset              nblocks;
    struct read_db_item     *next;
    struct read_db_item     *prev;
}read_db_item_t;

typedef struct ioreq_db{
    int                  updated_flag;
    MPI_Offset           nritems;
    struct write_db_item *witems;
    struct read_db_item  *ritems;
}ioreq_db_t;

typedef struct master_info{
    unsigned int         nrranks_completed;  /*Number of reader ranks that completed their read request. Counted only by master 0.*/
    unsigned int         nwranks_completed;  /*Number of writer ranks that completed their read requests
                                              (writer can also read the file). Counted only by master 0.*/
    unsigned int         nrranks_opened;     /*Number of reader ranks that opened the file*/
    unsigned int         nwranks_opened;     /*Number of writer ranks that opened the file*/
    struct ioreq_db      *iodb;
    int is_master_flag;  /*is this rank a master rank*/
    int nmasters;   /*Number of master nodes that hold data for request matching*/
    int *masters;   /*Ranks of master nodes on the writer's side*/
    unsigned int my_workgroup_sz;
} master_info_t;

typedef struct file_info_req_q{
    char filename[MAX_FILE_NAME];
    void *buf;  /*consists of
                  - filename [char[]]
                  - root reader rank to whom file info should be sent [int]
                  - number of readers that opened the file [int]*/
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




void add_ioreq(io_req_t **list, io_req_t *ioreq);
void delete_ioreqs(file_buffer_t *fbuf);
void progress_io_matching();
//void send_ioreq(file_buffer_t *fbuf, io_req_t *ioreq);
void send_ioreqs(file_buffer_t *fbuf);
void clean_iodb(ioreq_db_t *iodb);
int  match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag);
void match_ioreqs_all(int rw_flag);
int  init_req_match_masters(MPI_Comm comm, master_info_t *mst_info);
void init_iodb(file_buffer_t *fbuf);
void unpack_file_info(MPI_Offset bufsz, void *buf);
void send_file_info(file_buffer_t *fbuf, int reader_root);
#endif // dtf_REQ_MATCH_H_INCLUDED
