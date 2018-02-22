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
 
#define RDR_NOT_NOTIFIED   1
#define RDR_NOTIF_POSTED   2
#define RDR_NOTIFIED       3


typedef struct master_info{
    unsigned int         nread_completed;   /*Number of work groups that completed reading file*/
    struct ioreq_db      *iodb;
    int my_mst;     /*My master rank*/
    int nmasters;   /*Number of master nodes that hold data for request matching*/
    int *masters;   /*Ranks of master nodes on the writer's side*/
    unsigned int my_wg_sz;
    int *my_wg;
} master_info_t;

typedef struct dtf_var{
    int                     id;         /* varid assigned by pnetcdf*/
    MPI_Datatype            dtype;      /*Datatype of the variable*/
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
    struct io_req           *ioreqs;           /*Read or write I/O requests*/
    int                     max_dim;    /*Along what dimension the variable is broken into subblocks?
                                          Each master is responsible for metadata of a predefined subblock.
                                          */
    double                  checksum;   /*to check if what was written to the var is what was read*/
 }dtf_var_t;


typedef struct file_buffer{
  char                      file_path[MAX_FILE_NAME];    /* path of the file */
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
  int                       iomode;    /*Do normal File I/O or direct data transfer?*/
  int                       omode;     /*open mode (read/write/undefined)*/      
  int 						cur_transfer_epoch;      

  int                       root_writer;           /*MPI_COMM_WORLD rank of the rank who is a root in comm*/
  int                       root_reader;
  struct master_info        *cpl_mst_info;     /*Master info of the coupled component*/
  struct master_info        *my_mst_info;      /*Master info of this component*/
  struct io_req_log         *ioreq_log;        /*Used for debugging*/
  unsigned int              rreq_cnt;
  unsigned int              wreq_cnt;
  int                       done_matching_flag;     	/*Flag used to complete matching requests*/
  int                       sync_comp_flag;		/*Used for syncing two components*/
  int                       fready_notify_flag;   /*flag used to notify the reader that the file is ready for reading.
                                                    Possible values: 0 - the ps does not write to this file
                                                                     1 - the ps writes to this file. The file is not ready yet.
                                                                     2 - The file is ready, reader has been notified */
  int                       cpl_info_shared;                                          
  struct file_buffer        *next;
  struct file_buffer        *prev;
}file_buffer_t;

typedef struct fname_pattern{
    char fname[MAX_FILE_NAME];   /*File name pattern or the file name*/
    char **excl_fnames;         /*File name patterns to exclude*/
    int  nexcls;                 /*Number of file name patterns to exclude*/
    int  comp1;
    int  comp2;
    int  iomode;
    int  write_only;                
    int  ignore_io;				/*disable I/O for this file*/
    int  replay_io;            /*Should we record and replay during matching I/O pattern for this file?*/
    int  rdr_recorded;          /*if recorded, reader doesn't need to send its requests during matching as 
								the writer has already recorded its read pattern*/
    int  wrt_recorded;          /*If recorded, writer skips matching and replays data sending from previously 
                                 recorded I/O pattern*/
    struct io_pattern *io_pats;          /*Recorded pattern of what data this process sends to which reader process*/
    struct fname_pattern *next;
}fname_pattern_t;

void delete_file_buffer(file_buffer_t* buf);
file_buffer_t *create_file_buffer(fname_pattern_t *pat, const char* file_path);
fname_pattern_t* new_fname_pattern();
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid);
dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void add_var(file_buffer_t *fbuf, dtf_var_t *var);
int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start,const MPI_Offset *count );
fname_pattern_t *find_fname_pattern(const char *filename);

#endif
