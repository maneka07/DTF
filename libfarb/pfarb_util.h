#ifndef FARB_UTIL_H_INCLUDED
#define FARB_UTIL_H_INCLUDED

#include <stdlib.h>
#include <errno.h>
#include "pfarb_file_buffer.h"

#define FILE_READY_TAG  0
#define RECV_READY_TAG  1

#define HEADER_TAG  3
#define VARS_TAG    4
#define NODE_TAG    5 /*This tag should be always last!!*/

#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0

#define MAX_COMP_NAME 32
#define ASCIILINESZ   1024

/*
    DISTR_MODE_STATIC - configure data distribution through config file and
                        farb_set_distr_count(). Buffering happens on both sides:
                        reader's and writer's.

    DISTR_MODE_BUFFERED_REQ_MATCH - writer buffers all data, req matching
                        happens when the file is closed on writer side. Thanks to buffering,
                        it does not matter in what order the reader posts read requests and if this
                        order mismatches the write calls on the writer side. A deadlock should never happen.

    DISTR_MODE_NONBUFFERED_REQ_MATCH - writer does not buffer data, request matching happens inside wait().
                        Reader and writer are obliged to match the order of write and read requests, otherwise
                        a deadlock can happen. Reader must notify all writers that all its read requests have been
                        matched in order for the writer to return from a wait function.
*/
#define DISTR_MODE_STATIC                   0
#define DISTR_MODE_BUFFERED_REQ_MATCH       1
#define DISTR_MODE_NONBUFFERED_REQ_MATCH    2

//TODO: remove all asserts in the code and replace with proper error handling
typedef struct component{
    unsigned int id;
    //TODO remove +1
    char name[MAX_COMP_NAME+1];
    int connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    int master;       /*To whom to send read requests*/
    MPI_Comm intercomm;
}component_t;

typedef struct farb_config{
    int distr_mode;
    int my_master;  /* to whom I should send data about written data*/
    int nmasters;   /*Number of master nodes that hold data for request matching*/
    int *masters;   /*Ranks of master nodes on the writer's side*/
    unsigned int my_workgroup_sz;
}farb_config_t;



int get_read_flag(const char* filename);
int get_write_flag(const char* filename);
int get_io_mode(const char* filename);
void progress_io();
void notify_file_ready(const char* filename);
void close_file(const char* filename);
int file_buffer_ready(const char* filename);
void write_hdr(const char *filename, MPI_Offset hdr_sz, void *header);
MPI_Offset read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
int def_var(const char* filename, int varid, int ndims, MPI_Offset *shape);

int set_distr_count(const char* filename, int varid, int count[]);


MPI_Offset to_1d_index(int ndims, const MPI_Offset *shape, const MPI_Offset *coord);
#endif
