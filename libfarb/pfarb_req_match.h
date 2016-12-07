#ifndef PFARB_REQ_MATCH_H_INCLUDED
#define PFARB_REQ_MATCH_H_INCLUDED

#include <mpi.h>
#include "pfarb_file_buffer.h"

typedef struct contig_mem_chunk{
    MPI_Offset              offset;
    MPI_Offset              usrbuf_offset;
    MPI_Offset              data_sz;
    struct contig_mem_chunk *next;
}contig_mem_chunk_t;

typedef struct io_req{
    int                     id;
    int                     var_id;
    int                     completed;
    void                    *user_buf;
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
typedef struct write_chunk_rec{
    int                    rank;
    MPI_Offset             offset;
    MPI_Offset             data_sz;
    struct write_chunk_rec *next;
    struct write_chunk_rec *prev;
}write_chunk_rec_t;

typedef struct write_db_item{
    int                    var_id;
    struct write_chunk_rec *chunks;
    struct write_chunk_rec *last;
    struct write_db_item   *next;
}write_db_item_t;

typedef struct read_chunk_rec{
    int                     var_id;
    MPI_Offset              offset;
    MPI_Offset              data_sz;
    struct read_chunk_rec   *next;
    struct read_chunk_rec   *prev;
}read_chunk_rec_t;

typedef struct read_db_item{
    int                     rank;
    MPI_Offset              nchunks;
    struct read_chunk_rec   *chunks;
    struct read_db_item     *next;
    struct read_db_item     *prev;
}read_db_item_t;

typedef struct master_db{
    unsigned int         nranks_completed;
    unsigned int nmst_completed;
    MPI_Offset           nritems;
    struct write_db_item *witems;
    struct read_db_item  *ritems;
}master_db_t;

io_req_t *new_ioreq(int id,
                    int var_id,
                    int ndims,
                    MPI_Offset el_sz,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    void *buf,
                    int rw_flag,
                    int buffered);


void add_ioreq(io_req_t **list, io_req_t *ioreq);
void progress_io_matching();
void send_file_info(file_buffer_t *fbuf);
void send_ioreq(int ncid, io_req_t *ioreq, int rw_flag);
void clean_iodb(master_db_t *iodb);
int  match_ioreqs(file_buffer_t *fbuf);
void match_ioreqs_all(int rw_flag);

#endif // PFARB_REQ_MATCH_H_INCLUDED
