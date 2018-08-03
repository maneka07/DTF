#ifndef FILE_BUFFER_H_INCLUDED
#define FILE_BUFFER_H_INCLUDED

#include <mpi.h>
#include "rb_red_black_tree.h"
#include "dtf_var.h"

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

/*Structures related to I/O database that stores read and
 * writ I/O requests*/
typedef struct write_db_item{
    int                    ndims;
    void                   *dblocks; /*rb_red_blk_tree for ndims>0, block_t for ndims=0*/
    MPI_Offset             nblocks;
}write_db_item_t;

typedef struct read_dblock{
    int                     var_id;
    int                     ndims;
    MPI_Offset              *start;
    MPI_Offset              *count;
    struct read_dblock   *next;
    struct read_dblock   *prev;
}read_dblock_t;

typedef struct read_db_item{
    int                     global_rank;
    read_dblock_t           *dblocks;
//    read_dblock_t           *last_block;
    MPI_Offset              nblocks;
//    struct read_db_item     *next;
//    struct read_db_item     *prev;
}read_db_item_t;

typedef struct ioreq_db{
    int                  updated_flag;
    MPI_Offset           nritems;
    struct write_db_item **witems;
    struct read_db_item  **ritems;
}ioreq_db_t;

typedef struct master_info{
    unsigned int        	nread_completed;   /*Number of work groups that completed reading file*/
    struct ioreq_db     	*iodb;
    int 					my_mst;     /*My master rank*/
    int 					nmasters;   /*Number of master nodes that hold data for request matching*/
    int 					*masters;   /*Ranks of master nodes on the writer's side*/
    int                     comm_sz;    /*How many ranks opened the file*/
    unsigned int 			my_wg_sz;
} master_info_t;


typedef struct file_buffer{
  char                      file_path[MAX_FILE_NAME];    /* path of the file */
  int                       ncid;        /*handler that pnetcdf assigns to a file*/
  MPI_Comm                  comm;        /*MPI_Communicator used to open the file*/
  void                      *header;     /*buffer to store netcdf header*/
  MPI_Offset                hdr_sz;      /*size of the netcdf header*/
  dtf_var_t                 **vars;      
  int                       nvars;       /*Number of defined variables*/
  int                       writer_id;
  int                       reader_id;
  int                       is_ready;    /*Used to let the reader know that the file is either
										- received file header and var info from the writer (mode = DTF_IO_MODE_MEMORY)
										- finished being written (mode = DTF_IO_MODE_FILE)
										*/
  int                       is_transferring;	   /*Set to 1 when active transfer phase is happening*/
  int 						is_write_only;     /*If set, two components only write to this file.*/
  int                       iomode;            /*Do normal File I/O or direct data transfer?*/ 
  int                       is_defined;        /*is the file structure already defined?*/   //TODO use is_ready/is_transferring instead?? 
  int 						cur_transfer_epoch;      
  int                       root_writer;       /*MPI_COMM_WORLD rank of the rank who is a root in comm*/
  int                       root_reader;
  struct master_info        *cpl_mst_info;     /*Master info of the coupled component*/
  struct master_info        *my_mst_info;      /*Master info of this component*/

  unsigned long             rreq_cnt;
  unsigned long             wreq_cnt;
  int                       done_matching_flag;  /*Flag used to complete matching requests*/
  int                       done_multiple_flag;  /*Flag used to complete multiple transferes for the same file*/
  int                       fready_notify_flag;  /*flag used to notify the reader that the file is ready for reading.
                                                    Possible values: 0 - the ps does not write to this file
                                                                     1 - the ps writes to this file. The file is not ready yet.
                                                                     2 - The file is ready, reader has been notified */
  int 						session_cnt;         /*Every open/close is a session*/						
  int                       cpl_info_shared;    
  struct io_req_log         *ioreq_log;        /*Used for debugging*/
                                          
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
    int  num_sessions;          /*After how many open/close sessions we can delete given file buffer*/
    int  write_only;                
    int  replay_io;            /*Should we record and replay during matching I/O pattern for this file?*/
    int  rdr_recorded;          /*if recorded, reader doesn't need to send its requests during matching as 
								the writer has already recorded its read pattern*/
    int  wrt_recorded;          /*If recorded, writer skips matching and replays data sending from previously 
                                 recorded I/O pattern*/
    int  finfo_sz; 
    void *finfo;                /*File info (vars, header etc.). TO be used during I/O replaying*/
    struct io_pattern *io_pats;          /*Recorded pattern of what data this process sends to which reader process*/
    struct fname_pattern *next;
}fname_pattern_t;

file_buffer_t*		create_file_buffer(fname_pattern_t *pat, const char* file_path, MPI_Comm comm);
void 				delete_file_buffer(file_buffer_t* buf);
file_buffer_t* 		find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid);
void 				finalize_files();
void 				clean_iodb(ioreq_db_t *iodb, int nvars, int cpl_comm_sz);
void 				open_file(file_buffer_t *fbuf, MPI_Comm comm);
void 				close_file(file_buffer_t *fbuf);
void 				write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header);
MPI_Offset 			read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk);
void 				pack_file_info(file_buffer_t *fbuf, MPI_Offset *bufsz, void **buf);
void 				unpack_file_info(MPI_Offset bufsz, void *buf, file_buffer_t *fbf);
fname_pattern_t*	new_fname_pattern();
fname_pattern_t*	find_fname_pattern(const char *filename);



#endif
