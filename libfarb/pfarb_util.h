#ifndef FARB_UTIL_H_INCLUDED
#define FARB_UTIL_H_INCLUDED

#include <mpi.h>
#include "pfarb_file_buffer.h"


#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0

#define MAX_COMP_NAME 32
#define ASCIILINESZ   1024

//TODO if we decide not to implement a mixed writer-side-buffering + reader-side-req-matching
//then these tag definitions should be moved to corresponding .c files
#define FILE_READY_TAG      0       /*writer rank -> reader rank*/
#define RECV_READY_TAG      1       /*reader rank -> writer rank*/
#define HEADER_TAG          3       /*writer -> reader*/
#define VARS_TAG            4       /*writer -> reader*/
#define IO_WRITE_REQ_TAG    5       /*writer rank -> master rank*/
#define IO_READ_REQ_TAG     6       /*reader -> master*/
#define IO_DATA_REQ_TAG     7       /*master -> writer*/
#define IO_DATA_TAG         8       /*writer -> reader*/
#define IO_DONE_TAG         9       /*reader->master, master->writers*/
#define IO_CLOSE_FILE_TAG       10      /*reader->master, master->writers*/
#define NODE_TAG            11      /*This tag should be always last!!*/

/*
    DISTR_MODE_STATIC - configure data distribution through config file and
                        farb_set_distr_count(). Buffering happens on both sides:
                        reader's and writer's.

    DISTR_MODE_REQ_MATCH - //TODO.
                        Reader and writer are obliged to match the order of write and read requests, otherwise
                        a deadlock can happen. Reader must notify all writers that all its read requests have been
                        matched in order for the writer to return from a wait function.
*/
#define DISTR_MODE_STATIC                   0
#define DISTR_MODE_REQ_MATCH    1

//TODO: remove all asserts in the code and replace with proper error handling
typedef struct component{
    unsigned int    id;
    //TODO remove +1
    char            name[MAX_COMP_NAME+1];
    int             connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    MPI_Comm        comm;   /*intra or inter component communicator*/
}component_t;

typedef struct farb_config{
    int distr_mode;
    int my_master;  /*is this rank a master rank*/
    int nmasters;   /*Number of master nodes that hold data for request matching*/
    int *masters;   /*Ranks of master nodes on the writer's side*/
    unsigned int my_workgroup_sz;
    int buffered_req_match;    /*Should we buffer the data if request matching is enabled?*/
}farb_config_t;




struct file_buffer;

void progress_io();
void notify_file_ready(file_buffer_t *fbuf);
void close_file(file_buffer_t *fbuf);
int file_buffer_ready(const char* filename);
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header);
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Offset el_sz, MPI_Offset *shape);
int set_distr_count(file_buffer_t *fbuf, int varid, int count[]);
void open_file(file_buffer_t *fbuf);
MPI_Offset to_1d_index(int ndims, const MPI_Offset *shape, const MPI_Offset *coord);
MPI_Offset last_1d_index(int ndims, const MPI_Offset *shape);
#endif
