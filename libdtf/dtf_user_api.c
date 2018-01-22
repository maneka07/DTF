/*
* API intended to be used in user applications
*/
#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <errno.h>
#include <unistd.h>

#include "dtf.h"
#include "dtf_init_finalize.h"
#include "dtf_util.h"
#include "dtf_nbuf_io.h"
#include "dtf_req_match.h"

int lib_initialized=0;

file_info_req_q_t *gl_finfo_req_q = NULL;
file_info_t *gl_finfo_list = NULL;
dtf_msg_t *gl_msg_q = NULL;

struct file_buffer* gl_filebuf_list = NULL;        /*List of all file buffers*/
struct fname_pattern *gl_fname_ptrns = NULL;    /*Patterns for file name*/
struct component *gl_comps = NULL;                 /*List of components*/
int gl_my_comp_id;                          /*Id of this compoinent*/
int gl_ncomp;                               /*Number of components*/
int gl_verbose;
int gl_my_rank;                         /*For debug messages*/
int gl_scale;
struct dtf_config gl_conf;                 /*Framework settings*/
struct stats gl_stats;
char *gl_my_comp_name = NULL;
void* gl_msg_buf = NULL;

/**
  @brief	Function to initialize the library. Should be called from inside
            the application before any other call to the library. Should be called
            after the MPI is initialized.
  @param	filename        Name of the library configuration file.
  @param    module_name     Name of the module calling the init function.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int dtf_init(const char *filename, char *module_name)
{
  //  char* conf_filepath;
    int err, mpi_initialized;
    char* s;
    int verbose;

    if(lib_initialized)
        return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        fprintf(stderr, "DTF Error: dtf cannot be initialized before MPI is initialized. Aborting...\n");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &gl_my_rank);

    if(strlen(module_name)>MAX_COMP_NAME){
        fprintf(stderr, "DTF Error: module name %s too long\n", module_name);
        fflush(stderr);
        exit(1);
    }

    gl_stats.malloc_size = 0;
    gl_stats.data_msg_sz = 0;
    gl_stats.ndata_msg_sent = 0;
    gl_stats.accum_match_time = 0;
    gl_stats.ndb_match = 0;
    gl_stats.walltime = MPI_Wtime();
    gl_stats.accum_comm_time = 0;
    gl_stats.accum_hdr_time = 0;
    gl_stats.nprogress_call = 0;
    gl_stats.nioreqs = 0;
    gl_stats.nbl = 0;
    gl_stats.ngetputcall = 0;
    gl_stats.timer_accum = 0;
    gl_stats.timer_start = 0;
    gl_stats.accum_dbuff_sz = 0;
    gl_stats.accum_dbuff_time = 0;
    gl_stats.accum_rw_var = 0;
    gl_stats.accum_progr_time = 0;
    gl_stats.accum_do_matching_time = 0;
    gl_stats.nfiles = 0;
    gl_stats.idle_time = 0;
    gl_stats.idle_do_match_time = 0;
    gl_stats.master_time = 0;
    gl_stats.iodb_nioreqs = 0;
    gl_stats.parse_ioreq_time = 0;

    gl_my_comp_name = (char*)dtf_malloc(MAX_COMP_NAME);
    assert(gl_my_comp_name != NULL);
    strcpy(gl_my_comp_name, module_name);

    s = getenv("DTF_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_ERROR_LEVEL;
    else
        gl_verbose = atoi(s);

    //during init only root will print out stuff
    if(gl_my_rank != 0){
        verbose = gl_verbose;
        gl_verbose = VERBOSE_ERROR_LEVEL;
    }

    s = getenv("DTF_SCALE");
	if(s != NULL)
		gl_scale = atoi(s);
	else
		gl_scale = 0;
		
	DTF_DBG(VERBOSE_DBG_LEVEL, "Init DTF");

    s = getenv("DTF_DETECT_OVERLAP");
    if(s == NULL)
        gl_conf.detect_overlap_flag = 0;
    else
        gl_conf.detect_overlap_flag = atoi(s);

    s = getenv("DTF_DATA_MSG_SIZE_LIMIT");
    if(s == NULL)
        gl_conf.data_msg_size_limit = DTF_DATA_MSG_SIZE_LIMIT;
    else
        gl_conf.data_msg_size_limit = atoi(s) * 1024;
        
    
	s = getenv("DTF_IODB_RANGE");
    if(s == NULL)
        gl_conf.iodb_range = -1;
    else
        gl_conf.iodb_range = (MPI_Offset)atoi(s);
        
	DTF_DBG(VERBOSE_DBG_LEVEL, "Data message size limit set to %d", gl_conf.data_msg_size_limit);

    assert(gl_conf.data_msg_size_limit > 0);

    gl_msg_buf = NULL;
    gl_fname_ptrns = NULL;
    gl_filebuf_list = NULL;
    gl_finfo_req_q = NULL;
    gl_msg_q = NULL;

    /*Parse ini file and initialize components*/
    err = load_config(filename, module_name);
    if(err) goto panic_exit;

    /*Establish intercommunicators between components*/
    err = init_comp_comm();
    if(err) goto panic_exit;

    lib_initialized = 1;

    //enable print setting for other ranks again
    if(gl_my_rank != 0)
        gl_verbose = verbose;

    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF: Finished initializing");

    return 0;

panic_exit:

    dtf_finalize();
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    return 1;
}

/**
  @brief	Function to finalize the library. Should be called from inside
            the application before the MPI is finalized.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int dtf_finalize()
{
    int mpi_initialized, err;
    file_info_t *finfo;

    if(!lib_initialized) return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: dtf cannot be finalized after MPI is finalized. Aborting...");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }

//int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
//               MPI_Op op, int root, MPI_Comm comm)
    char *s = getenv("DTF_SCALE");
    if(s == NULL)
		DTF_DBG(VERBOSE_DBG_LEVEL,"time_stamp DTF: finalize");
	else 
		DTF_DBG(VERBOSE_ERROR_LEVEL,"time_stamp DTF: finalize");
    gl_stats.st_fin = MPI_Wtime() - gl_stats.walltime;
	
    /*Send any unsent file notifications
      and delete buf files*/
	while(gl_msg_q != NULL)
		progress_msg_queue();
		
    finalize_files();
    
    if(gl_msg_q != NULL)
		DTF_DBG(VERBOSE_DBG_LEVEL, "Finalize message queue");
	
    while(gl_msg_q != NULL)
		progress_msg_queue();

    assert(gl_finfo_req_q == NULL);
   
	finfo = gl_finfo_list;
	while(finfo != NULL){
		gl_finfo_list = gl_finfo_list->next;
		dtf_free(finfo, sizeof(file_info_t));
		finfo = gl_finfo_list;
	}
	
    finalize_comp_comm();
    print_stats();
    //destroy inrracomp communicator
    err = MPI_Comm_free(&gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);


    clean_config();

    if(gl_msg_buf != NULL)
        dtf_free(gl_msg_buf, gl_conf.data_msg_size_limit);

    //if(gl_stats.malloc_size != MAX_COMP_NAME )
      //  DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: DTF memory leak size: %lu", gl_stats.malloc_size - MAX_COMP_NAME);
    assert(gl_finfo_req_q == NULL);
    assert(gl_msg_q == NULL);

	if(gl_my_rank == 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF: finalized");
    dtf_free(gl_my_comp_name, MAX_COMP_NAME);
    lib_initialized = 0;
    fflush(stdout);
    fflush(stderr);
    return 0;
}

_EXTERN_C_ int dtf_match_io_v2(const char *filename, int ncid, int intracomp_io_flag, int it )
{
	char *s = getenv("DTF_IGNORE_ITER");
	if(it < 0)
		return dtf_match_io(filename, ncid, intracomp_io_flag);
		
	if(s != NULL){
		if(it > atoi(s)){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Match io call for iter %d", it);
			return dtf_match_io(filename, ncid, intracomp_io_flag);
		} else 
			DTF_DBG(VERBOSE_DBG_LEVEL, "Ignore match io call for iter %d", it);
	} else 
		return dtf_match_io(filename, ncid, intracomp_io_flag);
    return 0;
}

/*called by user to do explicit matching*/
/*
    User must specify either filename or ncid.
    intracomp_io_flag - if set to 1, matching of intracomponent io requests will be
    performed. This flag is intended for for situation when the writer component
    tries to read something from the file it is writing.
*/
//TODO rename to dtf_transfer()
_EXTERN_C_ int dtf_match_io(const char *filename, int ncid, int intracomp_io_flag )//, int match_all)
{
    file_buffer_t *fbuf;

    if(!lib_initialized) return 0;
    DTF_DBG(VERBOSE_DBG_LEVEL, "call match io for %s (ncid %d), intra flag %d", filename, ncid, intracomp_io_flag);
    if(intracomp_io_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: scale-letkf hack: skip intracomp matching");
        return 0;
    }

    fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
    if(fbuf == NULL){

        if( (filename != NULL) && (strlen(filename) == 0) )
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file (%s) with ncid %d is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        else
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        return 0;
    }
    if((fbuf->ioreq_log != NULL) && gl_conf.do_checksum){
		io_req_log_t *ioreq = fbuf->ioreq_log;
		while(ioreq != NULL){
			
			if(ioreq->dtype == MPI_DOUBLE || ioreq->dtype == MPI_FLOAT){
				double checksum = compute_checksum(ioreq->user_buf, ioreq->ndims, ioreq->count, ioreq->dtype);
				if(ioreq->rw_flag == DTF_READ)
					DTF_DBG(VERBOSE_ERROR_LEVEL, "read req %d, checksum %.4f", ioreq->id, checksum);
				else
					DTF_DBG(VERBOSE_ERROR_LEVEL, "write req %d, checksum %.4f", ioreq->id, checksum);
			}
			
			//delete
			fbuf->ioreq_log = ioreq->next;
			
			if(ioreq->rw_flag == DTF_READ)
				fbuf->rreq_cnt--;
			else
				fbuf->wreq_cnt--;
			
			if(gl_conf.buffered_req_match && (ioreq->rw_flag == DTF_WRITE)) dtf_free(ioreq->user_buf, ioreq->user_buf_sz);
			if(ioreq->start != NULL) dtf_free(ioreq->start, ioreq->ndims*sizeof(MPI_Offset));
			if(ioreq->count != NULL)dtf_free(ioreq->count, ioreq->ndims*sizeof(MPI_Offset));
			dtf_free(ioreq, sizeof(io_req_log_t));
			ioreq = fbuf->ioreq_log;
		}
	}
	 
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;

    if( intracomp_io_flag && (gl_my_comp_id != fbuf->writer_id) ){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_match_io: intracomp_io_flag(%d) can only be set for the writer component", intracomp_io_flag);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(fbuf->is_matching_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_io is called for file %s, but a matching process has already started before.", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

    dtf_tstart();
    match_ioreqs(fbuf, intracomp_io_flag);
    dtf_tend();
    return 0;
}

/*
 * This function is supposed to be called only by reader component to 
 * notify the writer component that it does not need the data. 
 * The function is added specifically for SCALE-LETKF because LETKF 
 * skips reading the history file when observation data for a given time 
 * slot is absent.*/
_EXTERN_C_ void dtf_skip_match(const char *filename, MPI_Comm comm)
{
	assert(0);
	//~ file_buffer_t *fbuf;
	//~ if(!lib_initialized) return;
	
	//~ if(filename == NULL || strlen(filename)==0){
		//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: undefined file name in dtf_skip_match call");
		//~ return;
	//~ }
	
	//~ DTF_DBG(VERBOSE_DBG_LEVEL, "Call skip match for %s", filename);
	
	//~ fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    //~ if(fbuf == NULL){
		//~ //If fbuf doesn't exist this file hasn't been opened/created. 
		//~ //However, we still need to check if I/O matching is enabled for 
		//~ //this file and, if it is, notify the other component that we skip.
		//~ fname_pattern_t *pat = gl_fname_ptrns;
		//~ while(pat != NULL){
			//~ if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls)){
				//~ DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
				//~ break;
			//~ }
			//~ pat = pat->next;
		//~ } 
		
		//~ if(pat != NULL){
			//~ if(gl_my_comp_id == pat->wrt){
				//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_skip_match for file %s should not be called by the writer component. Ignore. ", filename);
				//~ return;
			//~ }
			//~ assert(pat->rdr == gl_my_comp_id);
			
			//~ if(pat->iomode == DTF_IO_MODE_MEMORY){	
				//~ assert(comm != MPI_COMM_NULL);
				//~ skip_match(NULL, filename, comm, pat->wrt);
			//~ } 	
		//~ }
       
    //~ } else if(fbuf->iomode == DTF_IO_MODE_MEMORY)
		//~ skip_match(fbuf, filename, comm, -1);
	
}

_EXTERN_C_ void dtf_print(const char *str)
{
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_ERROR_LEVEL, "%s", str);
}


/*User controled timers*/
_EXTERN_C_ void dtf_time_start()
{
    if(!lib_initialized) return;
  
    //~ if(gl_stats.user_timer_start != 0)
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: user timer was started at %.3f and not finished.",
                //~ gl_stats.user_timer_start - gl_stats.walltime);

    gl_stats.user_timer_start = MPI_Wtime();
    DTF_DBG(VERBOSE_DBG_LEVEL, "user_time start");
    
    
    //~ char *s = getenv("DTF_SCALE");
	//~ if(s != NULL){
		//~ if(st_time_cnt==0){//strstr(fbuf->file_path, "hist.d")!=NULL)
			//~ gl_stats.st_mtch_hist = MPI_Wtime()-gl_stats.walltime;
		//~ } else if(st_time_cnt==1)//strstr(fbuf->file_path, "anal.d")!=NULL)
			//~ gl_stats.st_mtch_rest = MPI_Wtime()-gl_stats.walltime;
		//~ st_time_cnt++;
	//~ }
    
}
_EXTERN_C_ void dtf_time_end()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_stats.user_timer_start;
    gl_stats.user_timer_accum += tt;
    gl_stats.user_timer_start = 0;
    
	DTF_DBG(VERBOSE_DBG_LEVEL, "user_time end  %.6f", tt);
    
	//~ if(s != NULL){
		//~ if(end_time_cnt==0){//strstr(fbuf->file_path, "hist.d")!=NULL)
			//~ gl_stats.end_mtch_hist = MPI_Wtime()-gl_stats.walltime;
		//~ } else if(end_time_cnt==1)//strstr(fbuf->file_path, "anal.d")!=NULL)
			//~ gl_stats.end_mtch_rest = MPI_Wtime()-gl_stats.walltime;
		//~ end_time_cnt++;
	//~ }
    
 //   DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat: user time %.4f", tt);
}

/************************************************  Fortran Interfaces  *********************************************************/

_EXTERN_C_ void dtf_time_start_()
{
    dtf_time_start();
}
_EXTERN_C_ void dtf_time_end_()
{
    dtf_time_end();
}

void dtf_init_(const char *filename, char *module_name, int* ierr)
{
    *ierr = dtf_init(filename, module_name);
}

void dtf_finalize_(int* ierr)
{
    *ierr = dtf_finalize();
}

void dtf_match_io_(const char *filename, int *ncid, int *intracomp_io_flag, int *ierr)
{
    *ierr = dtf_match_io(filename, *ncid, *intracomp_io_flag);
}

void dtf_skip_match_(const char *filename, MPI_Fint *fcomm)
{
	return;
	dtf_skip_match(filename, MPI_Comm_f2c(*fcomm));
}	

void dtf_print_(const char *str)
{
    dtf_print(str);
}

void dtf_print_data_(int *varid, int *dtype, int *ndims, MPI_Offset* count, void* data)
{
    dtf_print_data(*varid, *dtype, *ndims, count, data);
}

