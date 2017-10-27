#ifndef FILE_BUFFER_H_INCLUDED
#define FILE_BUFFER_H_INCLUDED

#include <mpi.h>
#include "rb_red_black_tree.h"
//#include "dtf_req_match.h"

#define MAX_FILE_NAME 1024

/*
 0 (DTF_UNDEFINED) - the ps does not write to this file
 1 - the ps writes to this file. The file is not ready yet.
 2 - The file is ready, reader has been notified
 */
#define RDR_HASNT_OPENED   0
#define RDR_NOT_NOTIFIED   1
#define RDR_NOTIF_POSTED   2
#define RDR_NOTIFIED       3


typedef struct dtf_var{
    int                     id;         /* varid assigned by pnetcdf*/
    MPI_Datatype            dtype;      /*Datatype of the variable*/
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
    int                     max_dim;    /*Along what dimension the variable is broken into subblocks?
                                          Each master is responsible for metadata of a predefined subblock.
                                          */
    double                  checksum;   /*to check if what was written to the var is what was read*/
 }dtf_var_t;


struct master_db;
struct io_req;
struct master_struc;

typedef struct file_buffer{
  char                      file_path[MAX_FILE_NAME];    /* path of the file */
  char                      *slink_name;   
  int                       ncid;            /*handler that pnetcdf assigns to a file*/
  MPI_Comm                  comm;            /*MPI_Communicator used to open the file*/
  void                      *header;         /*buffer to store netcdf header*/
  MPI_Offset                hdr_sz;          /*size of the netcdf header*/
  dtf_var_t                 **vars;
  int                       nvars;            /*Number of defined variables*/
  int                       writer_id;
  int                       reader_id;
  int                       is_ready;           /*Used to let the reader know that the file is either
                                      - received from the writer (mode = DTF_IO_MODE_MEMORY)
                                      - finished being written (mode = DTF_IO_MODE_FILE)
                                    */
  int                       iomode;               /*Do normal File I/O or direct data transfer?*/
  /*Data distribution through request matching*/
  int                       root_writer;           /*MPI_COMM_WORLD rank of the rank who is a root in comm*/
  int                       root_reader;
  struct master_info       *mst_info;
  int                       nwriters;               /*Number of processes writing to the file*/
  //TODO move ioreqs inside var structure
  struct io_req             *ioreqs;           /*Read or write I/O requests*/
  int                       explicit_match;   /*0 - request matching is initiated from inside of pnetcdf;
                                                1 - request matching is initiated by the user*/
  unsigned int              rreq_cnt;
  unsigned int              wreq_cnt;
  int                       done_match_multiple_flag;   /*Set to 1 when reader notifies writer*/
  int                       done_matching_flag;     	/*Flag used to complete matching requests*/
  int                       done_match_confirm_flag;	/*Is set by the reader when the writer confirms that it finished matching.
														  Needed for multi-iterative programs to prevent the reader from sending
														  * ioreqs from the next iteration before writer finished with previous iteration*/
  int                       is_matching_flag;   /*Set to 1 when process starts matching. Reset to 0 when matching is done.*/
  int                       rdr_closed_flag;           /*Flag set to 1 when reader closes the file*/
  int                       fready_notify_flag;   /*flag used to notify the reader that the file is ready for reading.
                                                    Possible values: 0 - the ps does not write to this file
                                                                     1 - the ps writes to this file. The file is not ready yet.
                                                                     2 - The file is ready, reader has been notified */
  struct file_buffer        *next;             
  struct file_buffer        *prev;
}file_buffer_t;

typedef struct fname_pattern{
    char fname[MAX_FILE_NAME];   /*File name pattern or the file name*/
    char *slink_name;
    char **excl_fnames;         /*File name patterns to exclude*/
    int  nexcls;                 /*Number of file name patterns to exclude*/
    int  wrt;
    int  rdr;
    int  iomode;
    int  expl_mtch;
    struct fname_pattern *next;
}fname_pattern_t;

void delete_file_buffer(file_buffer_t* buf);
file_buffer_t *create_file_buffer(fname_pattern_t *pat, const char* file_path);
fname_pattern_t* new_fname_pattern();
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid);
dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void add_var(file_buffer_t *fbuf, dtf_var_t *var);
int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start,const MPI_Offset *count );
#endif
