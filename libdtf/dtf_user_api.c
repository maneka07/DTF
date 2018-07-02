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
#include "dtf_config.h"
#include "dtf_util.h"
#include "dtf_file_buffer.h"
#include "dtf_req_match.h"
#include "dtf_component.h"
#include "dtf_config.h"

int lib_initialized=0;

file_info_req_q_t*	gl_finfo_req_q 	= NULL;
file_info_t*		gl_finfo_list 	= NULL;
file_buffer_t*		gl_filebuf_list = NULL;			/*List of all file buffers*/
fname_pattern_t* 	gl_fname_ptrns 	= NULL;    		/*Patterns for file name*/
component_t* 		gl_comps 		= NULL;         /*List of components*/
int gl_my_comp_id = -1;                          	/*Id of this compoinent*/
int gl_ncomp = 0;                               	/*Number of components*/
int gl_verbose = 0;                         		/*For debug messages*/
int gl_scale = 0;
int gl_my_rank;
struct dtf_config 	gl_conf;                 		/*DTF settings*/
struct stats 		gl_stats;
char*	gl_my_comp_name = NULL;
void* 	gl_msg_buf = NULL;


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
    double t_start;
    
    if(lib_initialized)
        return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        fprintf(stderr, "DTF Error: dtf cannot be initialized before MPI is initialized. Aborting...\n");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }
    
    t_start = MPI_Wtime();
    
    MPI_Comm_rank(MPI_COMM_WORLD, &gl_my_rank);

    if(strlen(module_name)>MAX_COMP_NAME){
        fprintf(stderr, "DTF Error: module name %s too long\n", module_name);
        fflush(stderr);
        exit(1);
    }

    gl_stats.malloc_size = 0;
    gl_stats.data_msg_sz = 0;
    gl_stats.ndata_msg_sent = 0;
    gl_stats.transfer_time = 0;
    gl_stats.ndb_match = 0;
    gl_stats.walltime = MPI_Wtime();
    gl_stats.t_comm = 0;
    gl_stats.t_hdr = 0;
    gl_stats.nprogress_call = 0;
    gl_stats.nioreqs = 0;
    gl_stats.nbl = 0;
    gl_stats.ngetputcall = 0;
    gl_stats.timer_accum = 0;
    gl_stats.timer_start = 0;
    gl_stats.accum_dbuff_sz = 0;
    gl_stats.accum_dbuff_time = 0;
    gl_stats.t_rw_var = 0;
    gl_stats.t_progr_comm = 0;
    gl_stats.t_do_match = 0;
    gl_stats.nfiles = 0;
    gl_stats.idle_time = 0;
    gl_stats.idle_do_match_time = 0;
    gl_stats.master_time = 0;
    gl_stats.iodb_nioreqs = 0;
    gl_stats.parse_ioreq_time = 0;
    gl_stats.st_mtch_hist = 0;
    gl_stats.st_mtch_rest = 0;
    gl_stats.t_open_hist = 0;
    gl_stats.t_open_rest = 0;
    gl_stats.t_mtch_hist = 0;
    gl_stats.t_mtch_rest = 0;
	gl_stats.dtf_time = 0;
	gl_stats.t_idle = MPI_Wtime();

    gl_my_comp_name = (char*)dtf_malloc(MAX_COMP_NAME);
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
        
    
    s = getenv("DTF_IOREQ_FREQ");
    if(s == NULL)
        gl_conf.t_send_ioreqs_freq = 0.01;
    else
        gl_conf.t_send_ioreqs_freq = atof(s);
        
    assert(gl_conf.t_send_ioreqs_freq > 0);
    
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
	DTF_DBG(VERBOSE_ERROR_LEVEL, "Time to init DTF %.3f",  MPI_Wtime() - t_start);

    return 0;

panic_exit:
	DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: failed to initialize the framework. Aborting..");
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
    int mpi_initialized, err, comp;
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

	DTF_DBG(VERBOSE_ERROR_LEVEL,"time_stamp DTF: finalize");
	MPI_Barrier(gl_comps[gl_my_comp_id].comm);
    gl_stats.st_fin = MPI_Wtime() - gl_stats.walltime;
		
    /*Send any unsent file notifications
      //~ and delete buf files*/
    //~ if(gl_out_msg_q != NULL){
		//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "Send msg queue not empty:");
		//~ dtf_msg_t *msg = gl_out_msg_q;
		//~ while(msg != NULL){
			//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "%p (tag %d)", (void*)msg, msg->tag);
			//~ msg = msg->next;
		//~ }
		//~ while(gl_out_msg_q != NULL){
			//~ //progress_send_queue();
			//~ progress_comm();
		//~ }
	//~ }
	
	progress_send_queue();
		
    finalize_files();
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Notify other components that I finalized");
    for(comp = 0; comp < gl_ncomp; comp++){
		if(comp == gl_my_comp_id)
			continue;
		int ncomm_my, ncomm_cpl; 
		MPI_Comm_size(gl_comps[gl_my_comp_id].comm, &ncomm_my);
		MPI_Comm_remote_size(gl_comps[comp].comm, &ncomm_cpl);
		if(gl_my_rank < ncomm_cpl){
			dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, COMP_FINALIZED_TAG);
			err = MPI_Isend(NULL, 0, MPI_INT, gl_my_rank, COMP_FINALIZED_TAG, gl_comps[comp].comm, &(msg->req));
			CHECK_MPI(err);
			ENQUEUE_ITEM(msg, gl_comps[comp].out_msg_q);
		}
		
		if( (ncomm_cpl > ncomm_my) && (gl_my_rank == 0)){
			int i;
			for(i = gl_my_rank; i < ncomm_cpl; i++){
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, COMP_FINALIZED_TAG);
				err = MPI_Isend(NULL, 0, MPI_INT, gl_my_rank, COMP_FINALIZED_TAG, gl_comps[comp].comm, &(msg->req));
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_comps[comp].out_msg_q);
			}
		} 
	}
    
    for(comp = 0; comp < gl_ncomp; comp++){
		if(gl_comps[comp].out_msg_q == NULL)
			continue;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Finalize message queue for comp %s", gl_comps[comp].name);
		while(gl_comps[comp].out_msg_q != NULL)
			progress_comm();	
	}
		
    
    assert(gl_finfo_req_q == NULL);
   
	finfo = gl_finfo_list;
	while(finfo != NULL){
		gl_finfo_list = gl_finfo_list->next;
		dtf_free(finfo, sizeof(file_info_t));
		finfo = gl_finfo_list;
	}
	
	gl_comps[gl_my_comp_id].finalized = 1;
    
    finalize_comp_comm();
    print_stats();
    //destroy inrracomp communicator
    err = MPI_Comm_free(&gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    
	for(comp = 0; comp < gl_ncomp; comp++)
		assert(gl_comps[comp].out_msg_q == NULL);

    clean_config();

    if(gl_msg_buf != NULL)
        dtf_free(gl_msg_buf, gl_conf.data_msg_size_limit);

  //  if(gl_stats.malloc_size != MAX_COMP_NAME )
    //  DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: DTF memory leak size: %lu", gl_stats.malloc_size - MAX_COMP_NAME);
    assert(gl_finfo_req_q == NULL);
    
	DTF_DBG(VERBOSE_DBG_LEVEL,"DTF: finalized");
    dtf_free(gl_my_comp_name, MAX_COMP_NAME);
    lib_initialized = 0;
    fflush(stdout);
    fflush(stderr);
    return 0;
}

_EXTERN_C_ int dtf_transfer_v2(const char *filename, int ncid, int it )
{
	if(!lib_initialized) return 0;
	char *s = getenv("DTF_IGNORE_ITER");
	
	if(it < 0)
		return dtf_transfer(filename, ncid);
		
	if(s != NULL){
		if(it > atoi(s)){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Match io call for iter %d", it);
			return dtf_transfer(filename, ncid);
		} else 
			DTF_DBG(VERBOSE_DBG_LEVEL, "Ignore match io call for iter %d", it);
	} else 
		return dtf_transfer(filename, ncid);
    return 0;
}

/*
 * Start non-blocking data transfer. Processes will try to progress
 * with data transfer as much as they can but if there is no more work 
 * to do at the moment, the process will simply return. 
 * The user needs to eventually call dtf_transfer_complete or 
 * dtf_transfer_complete_all to ensure the completion of the data tranfer.
 * 
 * NOTE: The user cannot call dtf_transfer_start or dtf_transfer for the 
 * same file if there is already an active data transfer
 * */
//_EXTERN_C_ int dtf_transfer_start(const char *filename)
//{
	
//}

/*
 * Process will block until all active data transfers are completed.
 * */
_EXTERN_C_ int dtf_transfer_all_files()
{
	double t_start;
	if(!lib_initialized) return 0;
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Start transfer_complete_all");
	t_start = MPI_Wtime();
	match_ioreqs_all_files();
	gl_stats.transfer_time += MPI_Wtime() - t_start;
	gl_stats.dtf_time += MPI_Wtime() - t_start;
	DTF_DBG(VERBOSE_DBG_LEVEL, "End transfer_complete_all");
	DTF_DBG(VERBOSE_DBG_LEVEL, "Time complall %.3f",  MPI_Wtime() - t_start);
	return 0;
}


/*called by user to do explicit matching*/
/*
    User must specify either filename or ncid.
*/

_EXTERN_C_ int dtf_transfer(const char *filename, int ncid)
{
    file_buffer_t *fbuf;
	double t_start = MPI_Wtime();
	
    if(!lib_initialized) return 0;
    DTF_DBG(VERBOSE_DBG_LEVEL, "call match io for %s (ncid %d)", filename, ncid);
  
    fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
    if(fbuf == NULL){

        if( (filename != NULL) && (strlen(filename) == 0) )
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file (%s) with ncid %d is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        else
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        return 0;
    }
	 
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    
	if(!fbuf->is_transferring) return 0;
    
    match_ioreqs(fbuf);
     
	fname_pattern_t *pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions)
		delete_file_buffer(fbuf);
    
    gl_stats.transfer_time += MPI_Wtime() - t_start;
    gl_stats.dtf_time += MPI_Wtime() - t_start;
    return 0;
}

/*Supposed to be called by the writer process.
  Used to match against several dtf_match_io functions on the reader side*/
_EXTERN_C_ void dtf_transfer_multiple(const char* filename, int ncid)
{
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s ncid %d is not treated by DTF( \
                not in configuration file). Explicit matching ignored.", filename, ncid);
        return;
    }
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(!fbuf->is_transferring) return;
    
	double t_start = MPI_Wtime();
    if( gl_my_comp_id != fbuf->writer_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_transfer_multiple can only be called by writer component.");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching multiple");
   
    fbuf->done_multiple_flag = 0;
    while(!fbuf->done_multiple_flag){
        match_ioreqs(fbuf);
        
        /*If the other component has started finalizing, then just complete this transfer*/
        if( ((fbuf->writer_id == gl_my_comp_id) && gl_comps[fbuf->reader_id].finalized)  ||
			((fbuf->reader_id == gl_my_comp_id) && gl_comps[fbuf->writer_id].finalized) ){
				fbuf->done_multiple_flag = 1;
		}
	}
    //reset
    fbuf->done_multiple_flag = 0;
    
	fname_pattern_t *pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions)
		delete_file_buffer(fbuf);
    
    gl_stats.transfer_time += MPI_Wtime() - t_start;
    gl_stats.dtf_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finish matching multiple");
    return;
}

/*Used by reader to notify writer that it can complete dtf_match_multiple*/
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid)
{
	fname_pattern_t *pat;
    double t_start = MPI_Wtime();
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). matching ignored.", filename, ncid);
        return;
    }

    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    
    if(fbuf->is_write_only) return;
    
    if(fbuf->reader_id != gl_my_comp_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_complete_multiple for file %s can only be called by reader", filename);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    
    MPI_Barrier(fbuf->comm); //TODO remove barrier?
    
    if( ((fbuf->writer_id == gl_my_comp_id) && !gl_comps[fbuf->reader_id].finalized)  ||
		((fbuf->reader_id == gl_my_comp_id) && !gl_comps[fbuf->writer_id].finalized) )			
		notify_complete_multiple(fbuf);
    
    pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions && !fbuf->is_transferring)
		delete_file_buffer(fbuf);
    
    gl_stats.transfer_time += MPI_Wtime() - t_start;
    gl_stats.dtf_time += MPI_Wtime() - t_start;
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
  
    gl_stats.user_timer_start = MPI_Wtime();
    DTF_DBG(VERBOSE_DBG_LEVEL, "user_time start");
    
    
}

_EXTERN_C_ void dtf_time_end()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_stats.user_timer_start;
    gl_stats.user_timer_accum += tt;
    gl_stats.user_timer_start = 0;
    
	DTF_DBG(VERBOSE_DBG_LEVEL, "user_time end  %.6f", tt);
    
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

void dtf_transfer_(const char *filename, int *ncid, int *ierr)
{
    *ierr = dtf_transfer(filename, *ncid);
}

void dtf_transfer_v2_(const char *filename, int *ncid, int *it, int *ierr )
{
	*ierr = dtf_transfer_v2(filename, *ncid, *it);
}

void dtf_transfer_all_files_()
{
	dtf_transfer_all_files();
}

void dtf_transfer_multiple_(const char *filename, int *ncid)
{
    dtf_transfer_multiple(filename, *ncid);
}

void dtf_complete_multiple_(const char *filename, int *ncid)
{
    dtf_complete_multiple(filename, *ncid);
}

void dtf_print_(const char *str)
{
    dtf_print(str);
}

void dtf_print_data_(int *varid, int *dtype, int *ndims, MPI_Offset* count, void* data)
{
    dtf_print_data(*varid, *dtype, *ndims, count, data);
}

