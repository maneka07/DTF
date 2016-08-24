#ifndef FARB_UTIL_H_INCLUDED
#define FARB_UTIL_H_INCLUDED

#include <stdlib.h>
#include <errno.h>

#include "pfarb_common.h"
#include "pfarb_file_buffer.h"

#define FILE_READY_TAG  0
#define RECV_READY_TAG  1

#define HEADER_TAG  3
#define VARS_TAG    4
#define NODE_TAG    5 /*This tag should be always last!!*/

#define CONNECT_MODE_UNDEFINED -1
#define CONNECT_MODE_SERVER     1
#define CONNECT_MODE_CLIENT     0


//#include "farb_file_buffer.h"

//TODO: remove all asserts in the code and replace with proper error handling
//TODO file can have either only reader or writer
typedef struct component{
    unsigned int id;
    //TODO remove +1
    char name[MAX_COMP_NAME+1];
    int connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    MPI_Comm intercomm;
}component_t;

typedef struct farb_settings{
    MPI_Offset node_sz;         /*Size of a node of memory buffer*/
    int msg_sz;             /*MPI message size for transfering the data. Should be a divisor for node_sz.*/
}farb_settings_t;

int load_config(const char *ini_name, const char *service_name);
void clean_config();

int init_comp_comm();
void finalize_comp_comm();

int get_read_flag(const char* filename);
int get_write_flag(const char* filename);
int get_io_mode(const char* filename);
void progress_io();
void notify_file_ready(const char* filename);
void notify_recv_ready(const char* filename);
void close_file(const char* filename);
int file_buffer_ready(const char* filename);
void write_hdr(const char *filename, MPI_Offset hdr_sz, void *header);
void pack_vars(int var_cnt, farb_var_t *vars, int *buf_sz, void **buf);
void unpack_vars(int buf_sz, void *buf, farb_var_t **vars, int *var_cnt);

MPI_Offset read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
int def_var(const char* filename, int varid, int ndims, MPI_Offset *shape);
MPI_Offset read_write_var(const char *filename, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, MPI_Datatype dtype, void *buf, int rw_flag);
#endif
