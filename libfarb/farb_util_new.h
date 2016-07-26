#ifndef FARB_UTIL_H_INCLUDED
#define FARB_UTIL_H_INCLUDED

#include <stdlib.h>
#include <errno.h>

#include "farb_common.h"
#include "farb_file_buffer.h"

#define FILE_READY_TAG  0
#define RECV_READY_TAG  1
#define DATA_TAG    2

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
    size_t node_sz;         /*Size of a node of memory buffer*/
    int msg_sz;             /*MPI message size for transfering the data. Should be a divisor for node_sz.*/
}farb_settings_t;

typedef struct msg_ready_notif{
    char filename[MAX_FILE_NAME];
    unsigned int file_sz;
}msg_ready_notif_t;


MPI_Datatype msg_ready_datatype;

int load_config(const char *ini_name, const char *service_name);
//void print_config();
void clean_config();

int init_comp_comm();
void finalize_comp_comm();

size_t mem_write(const char* filename, off_t const offset, const size_t data_sz, void *const data);
size_t mem_read(const char* filename, off_t const offset, const size_t data_sz, void *const data);
int get_read_flag(const char* filename);
int get_write_flag(const char* filename);
int get_io_mode(const char* filename);
void progress_io();
void notify_file_ready(const char* filename);
void notify_recv_ready(const char* filename);
void close_file(const char* filename);
int file_buffer_ready(const char* filename);
#endif
