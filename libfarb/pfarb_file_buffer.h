#ifndef FILE_BUFFER_H_INCLUDED
#define FILE_BUFFER_H_INCLUDED

#include "pfarb_common.h"
#include "pfarb_buffer_node.h"
#include <string.h>

typedef struct farb_var{
    int                     id;         /* varid assigned by pnetcdf*/
 //   int                     xsz;        /* byte size of 1 array element */
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
//    int                     *dimids;    /* array of dimension IDs */
 //   MPI_Datatype            type;       /* variable's data type */
  //  MPI_Offset              begin;      /* starting file offset of this variable */
 //   MPI_Offset              last_byte_offt;    /*offset of the last byte written by this ps*/

    buffer_node_t           *nodes;      /*head buffer node that stores the data*/
    MPI_Offset              node_cnt;   /*number of allocated buffer nodes*/
 //    MPI_Offset              data_sz;    /*total amount of data of this variable that this ps wrote*/

    struct farb_var *next;
}farb_var_t;

typedef struct file_buffer{
  char file_path[MAX_FILE_NAME];    /* path of the file */
  char alias_name[MAX_FILE_NAME];	/* alias name for the file */

  void          *header;            /*buffer to store netcdf header*/
  MPI_Offset    hdr_sz;             /*size of the netcdf header*/
  struct farb_var *vars;              /*Variables in the file*/
  int           var_cnt;            /*Number of defined variables*/
  int           version;            /*To keep track in case the component generates
                                    several output files with the same name pattern*/
  int           writer_id;          /*Only one component can write to a file*/
  int           *reader_ids;        /*Multiple components can read from a file*/
  int           nreaders;           /*Number of components that read from the file*/
  int           is_ready;           /*Used to let the reader know that the file is either
                                      - received from the writer (mode = FARB_IO_MODE_MEMORY)
                                      - finished being written (mode = FARB_IO_MODE_FILE)
                                    */
  int           write_flag;         /*0 - this component does not write to memory for this file, otherwise - 1*/
  int           read_flag;          /*0 - this component does not read from memory for this file, otherwise - 1*/
  int           transfered;         /*Incremented every time  the writer either
                                      - sends the file to one of the readers(mode = FARB_IO_MODE_MEMORY)
                                      - notifies a reader that it finished writing the file (mode = FARB_IO_MODE_FILE)*/
  int           mode;               /*Do normal File I/O or direct data transfer?*/

  struct file_buffer *next;         /* pointer to the next record */

}file_buffer_t;

void add_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);
void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);

file_buffer_t* new_file_buffer();
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path);
farb_var_t* find_var(farb_var_t* varlist, int varid);
farb_var_t* new_var(int varid, int ndims, MPI_Offset *shape);
void add_var(farb_var_t **vars, farb_var_t *var);
#endif
