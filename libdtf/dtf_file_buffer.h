#ifndef FILE_BUFFER_H_INCLUDED
#define FILE_BUFFER_H_INCLUDED

#include "dtf_buffer_node.h"
//#include "dtf_req_match.h"

#define MAX_FILE_NAME 1024

#define DISTR_RULE_P2P      0  /*default*/
#define DISTR_RULE_RANGE    1

#define DISTR_PATTERN_SCATTER   0
#define DISTR_PATTERN_ALL       1

/*
 0 (DTF_UNDEFINED) - the ps does not write to this file
 1 - the ps writes to this file. The file is not ready yet.
 2 - The file is ready, reader has been notified
 */
 #define RDR_NOT_NOTIFIED   1
 #define RDR_NOTIFIED       2


typedef struct dtf_var{
    int                     id;         /* varid assigned by pnetcdf*/
    //MPI_Offset              el_sz;      /* byte size of 1 array element */
    MPI_Datatype            dtype;      /*Datatype of the variable*/
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
    buffer_node_t           *nodes;     /*head buffer node that stores the data*/
    MPI_Offset              node_cnt;   /*number of allocated buffer nodes*/
    MPI_Offset *first_coord;            /*Coordinate of the first element in this node. Needed for scatter distribution. Needed only by writer.*/
/*Static data distribution related stuff*/
    MPI_Offset *distr_count;                 /*Number of elements in each dimension to distribute*/
    struct dtf_var *next;
}dtf_var_t;


struct master_db;
struct io_req;
struct master_struc;

typedef struct file_buffer{
  char                      file_path[MAX_FILE_NAME];    /* path of the file */
  char                      alias_name[MAX_FILE_NAME];	/* alias name for the file */
  int                       ncid;                        /*handler that pnetcdf assigns to a file*/
  MPI_Comm                  comm;               /*MPI_Communicator used to open the file*/
  void                      *header;            /*buffer to store netcdf header*/
  MPI_Offset                hdr_sz;             /*size of the netcdf header*/
  struct dtf_var           *vars;              /*Variables in the file*/
  int                       var_cnt;            /*Number of defined variables*/
  int                       version;            /*To keep track in case the component generates
                                    several output files with the same name pattern*/
  int                       writer_id;
  int                       reader_id;
  int                       is_ready;           /*Used to let the reader know that the file is either
                                      - received from the writer (mode = DTF_IO_MODE_MEMORY)
                                      - finished being written (mode = DTF_IO_MODE_FILE)
                                    */
  int                       iomode;               /*Do normal File I/O or direct data transfer?*/
  /*Static data distribution related stuff*/
  int                       distr_rule;                   /*Rule for distribution from writer to readers(range or list of ranks)*/
  int                       distr_pattern;                /*Scatter portions of data or send all data*/
  int                       distr_range;                  /*Range for distribution*/
  int                       distr_nranks;                 /*writer: how many ranks I distribute to; reader: how many ranks I receive from*/
  int                       *distr_ranks;                 /*writer: ranks I distribute to; reader: ranks I receive from*/
  int                       distr_ndone;                  /*number of completed distributions*/
  /*Data distribution through request matching*/
  int                       root_writer;           /*MPI_COMM_WORLD rank of the rank who is a root in comm*/
  int                       root_reader;
  struct master_info       *mst_info;
  int                       nwriters;               /*Number of processes writing to the file*/
  struct io_req             *ioreqs;           /*Read or write I/O requests*/
  int                       explicit_match;   /*0 - request matching is initiated from inside of pnetcdf;
                                                1 - request matching is initiated by the user*/
  //int                       ioreq_cnt;         /*Request counter, will be used to assign a unique id to io requests.*/
  unsigned int              rreq_cnt;
  unsigned int              wreq_cnt;
  int                       done_matching_flag;     /*Flag used to complete matching requests*/
  int                       is_matching_flag;   /*Set to 1 when process starts matching. Reset to 0 when matching is done.*/
  int                       rdr_closed_flag;           /*Flag set to 1 when reader closes the file*/
  int                       fready_notify_flag;   /*flag used to notify the reader that the file is ready for reading.
                                                    Possible values: 0 - the ps does not write to this file
                                                                     1 - the ps writes to this file. The file is not ready yet.
                                                                     2 - The file is ready, reader has been notified */
  struct file_buffer        *next;             /* pointer to the next record */

}file_buffer_t;

void add_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);
void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf);

file_buffer_t* new_file_buffer();
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid);
dtf_var_t* find_var(dtf_var_t* varlist, int varid);
dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void add_var(dtf_var_t **vars, dtf_var_t *var);
int has_unlim_dim(dtf_var_t *var);
int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start,const MPI_Offset *count );
#endif
