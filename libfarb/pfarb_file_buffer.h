#ifndef FILE_BUFFER_H_INCLUDED
#define FILE_BUFFER_H_INCLUDED

#include "farb_common.h"
#include "farb_buffer_node.h"
#include <string.h>

typedef struct file_buffer{
  char file_path[MAX_FILE_NAME];         /* path of the file */
  char alias_name[MAX_FILE_NAME];		/* alias name for the file */

  buffer_node_t* nodes_tbl;            /*for faster access keep an array of pointers to nodes*/
  unsigned int node_cnt;
  size_t data_sz;                        /*size of the data stored in the buffer*/

  /*Fields for the pnetcdf implementation*/
  int version;                  /*To keep track in case the component generates
                                several output files with the same name pattern*/
  int writer_id;                /*Only one component can write to a file*/
  int* reader_ids;               /*Multiple components can read from a file*/
  int  nreaders;                 /*Number of components that read from the file*/
  int is_ready;                   /*Used to let the reader know that the file is either
                                    - received from the writer (mode = IO_MODE_MEMORY)
                                    - finished being written (mode = IO_MODE_FILE)
                                    */
  int  write_flag;               /*0 - this component does not write to memory for this file, otherwise - 1*/
  int  read_flag;                 /*0 - this component does not read from memory for this file, otherwise - 1*/
  int  transfered;                /*Incremented every time  the writer either
                                    - sends the file to one of the readers(mode = IO_MODE_MEMORY)
                                    - notifies a reader that it finished writing the file (mode = IO_MODE_FILE)*/
  int mode;                       /*Do normal File I/O or direct data transfer?*/

  struct file_buffer *next;     /* pointer to the next record */

  /*TODO: These will be deleted in the future*/
  buffer_node_t* nodes;                 /*list of nodes that hold the data*/
  buffer_node_t* last_node;
  char *buffer_pointer;         /* pointer of the buffer */
  struct buffer_node *buffer_list; /* buffer node list */
  int buffer_size;              /* size of the buffer */
  char service_name[1024];		/* service name */
  char direction[8];			/* direction of send & recv */
  int client_server_flag;
  int files_cnt;

}file_buffer_t;

void add_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);
void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);
file_buffer_t* new_file_buffer();
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path);
#endif
