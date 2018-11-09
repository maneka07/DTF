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
int gl_verbose = 0;                         		/*For debug messages*/
struct dtf_proc gl_proc;

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
    
    if(strlen(module_name)>MAX_COMP_NAME){
        fprintf(stderr, "DTF Error: module name %s too long\n", module_name);
        fflush(stderr);
        exit(1);
    }
    
    strcpy(_comp_name, module_name);
    
    gl_proc.stats_info.num_fpats = 0;
    gl_proc.stats_info.malloc_size = 0;
    gl_proc.stats_info.data_msg_sz = 0;
    gl_proc.stats_info.ndata_msg_sent = 0;
    gl_proc.stats_info.transfer_time = 0;
    gl_proc.walltime = MPI_Wtime();
    gl_proc.stats_info.t_comm = 0;
    gl_proc.stats_info.nioreqs = 0;
    gl_proc.stats_info.timer_accum = 0;
    gl_proc.stats_info.timer_start = 0;
    gl_proc.stats_info.accum_dbuff_sz = 0;
    gl_proc.stats_info.accum_dbuff_time = 0;
    gl_proc.stats_info.t_progr_comm = 0;
    gl_proc.stats_info.t_do_match = 0;
    gl_proc.stats_info.nfiles = 0;
    gl_proc.stats_info.idle_time = 0;
    gl_proc.stats_info.idle_do_match_time = 0;
    gl_proc.stats_info.master_time = 0;
    gl_proc.stats_info.iodb_nioreqs = 0;
	gl_proc.stats_info.dtf_time = 0;
	gl_proc.stats_info.t_idle = MPI_Wtime();
	gl_proc.msgbuf = NULL;
    gl_proc.fname_ptrns = NULL;
    gl_proc.filebuf_list = NULL;
	gl_proc.comps 		= NULL;       
	gl_proc.my_comp = DTF_UNDEFINED; 
	gl_proc.ncomps = 0;                               
	gl_proc.conf.log_ioreqs = 0;
    gl_proc.conf.buffer_data = 0;
    gl_proc.conf.do_checksum = 0;
    gl_proc.conf.iodb_build_mode = IODB_BUILD_BLOCK; // default     
    
    s = getenv("DTF_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_ERROR_LEVEL;
    else
        gl_verbose = atoi(s);

	MPI_Comm_rank(MPI_COMM_WORLD, &(gl_proc.myrank));

    //during init only root will print out stuff
    if(gl_proc.myrank != 0){
        verbose = gl_verbose;
        gl_verbose = VERBOSE_ERROR_LEVEL;
    }
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Init DTF");
		
    s = getenv("DTF_DATA_MSG_SIZE_LIMIT");
    if(s == NULL)
        gl_proc.conf.data_msg_size_limit = DTF_DATA_MSG_SIZE_LIMIT;
    else
        gl_proc.conf.data_msg_size_limit = atoi(s) * 1024;      
    
	s = getenv("DTF_USE_MSG_BUFFER");
    if(s == NULL)
		gl_proc.conf.use_msg_buffer = 0;
	else 
		gl_proc.conf.use_msg_buffer = atoi(s);
    
	s = getenv("DTF_VAR_BLOCK_RANGE");
    if(s == NULL)
        gl_proc.conf.iodb_range = 0;
    else
        gl_proc.conf.iodb_range = (MPI_Offset)atoi(s);
    
    s = getenv("DTF_TIMEOUT");
	if(s == NULL)
		gl_proc.conf.timeout = 180;
	else
		gl_proc.conf.timeout = atoi(s);
    
    s = getenv("DTF_IGNORE_IDLE");
    if(s == NULL)
        gl_proc.conf.ignore_idle = 0;
    else
        gl_proc.conf.ignore_idle = atoi(s);
        
	DTF_DBG(VERBOSE_DBG_LEVEL, "Data message size limit set to %d", gl_proc.conf.data_msg_size_limit);

    assert(gl_proc.conf.data_msg_size_limit > 0);

    
	/*Parse ini file and initialize components*/
	err = load_config(filename, module_name);
	if(err) goto panic_exit;

    /*find out if we run all components with single mpirun/mpiexec or 
     * separate command for each component. */
    
	{
		//Temporarily turn off profiling if it's running
		MPI_Pcontrol(0);
		
		void *val;
		int flag, vval, maxval=-1;

		err = MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_APPNUM, &val, &flag );
		CHECK_MPI(err);
		if(flag){
			vval = *(int*)val; 
			err = MPI_Allreduce(&vval, &maxval, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
			CHECK_MPI(err);
			DTF_DBG(VERBOSE_DBG_LEVEL, "APP_NUM is %d", maxval);
			gl_proc.conf.single_mpirun_mode =  maxval > 0 ? 1 : 0; 
		}else {
			
			DTF_DBG(VERBOSE_DBG_LEVEL, "APP_NUM is not set!!");
			char *s = getenv("MPMD_COMP");
			gl_proc.conf.single_mpirun_mode = s != NULL ? 1: 0;
		    
		}
		MPI_Pcontrol(1);
		
		if(gl_proc.conf.single_mpirun_mode)
			DTF_DBG(VERBOSE_DBG_LEVEL, "Consider single mpirun execution mode");
		else 
			DTF_DBG(VERBOSE_DBG_LEVEL, "Consider separate mpirun for each executable execution mode");
	}
	
        
    /*Establish intercommunicators between components*/
    err = init_comp_comm();
    if(err) assert(0); //goto panic_exit;

    create_tmp_file();
   
    lib_initialized = 1;

    //enable print setting for other ranks again
    if(gl_proc.myrank != 0)
        gl_verbose = verbose;

	if(gl_proc.myrank==0)
		DTF_DBG(VERBOSE_DBG_LEVEL, "Time to init DTF %.3f",  MPI_Wtime() - t_start);

    return 0;

panic_exit:
	DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Error: failed to initialize the framework. Aborting..");
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

    if(!lib_initialized) return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: dtf cannot be finalized after MPI is finalized. Aborting...");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }

	DTF_DBG(VERBOSE_DBG_LEVEL,"time_stamp DTF: finalize");
	MPI_Barrier(gl_proc.comps[gl_proc.my_comp].comm);
	
	progress_send_queue();
		
    finalize_files();
    
    delete_tmp_file();
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Notify other components that I finalized");
    if(gl_proc.conf.single_mpirun_mode){
		/*Send finalize message to ranks of other components before and 
		 * after my component
		 *       ----------------------------
		 * 		 | others   | Me   |  others|
		 * 		 ----------------------------
		 * */
		 //found out how many ranks before and after ranks of my component
		 int my_comm_sz, glob_comm_sz;
		 int me_comp, me_glob;
		 int num_msg, i, start_rank;
		 
		for(comp = 0; comp < gl_proc.ncomps; comp++){
			if(comp == gl_proc.my_comp)
				continue;
		 
			MPI_Comm_rank(gl_proc.comps[gl_proc.my_comp].comm, &me_comp);
			MPI_Comm_rank(gl_proc.comps[comp].comm, &me_glob);
			MPI_Comm_size(gl_proc.comps[gl_proc.my_comp].comm, &my_comm_sz);
			MPI_Comm_size(gl_proc.comps[comp].comm, &glob_comm_sz);
			
			num_msg = (glob_comm_sz - my_comm_sz)/my_comm_sz;
			if(num_msg == 0 && (me_comp < glob_comm_sz - my_comm_sz))
				num_msg = 1;
			if(glob_comm_sz - my_comm_sz > my_comm_sz && 
			   (glob_comm_sz - my_comm_sz)%my_comm_sz != 0 && 
			    me_comp == my_comm_sz-1)
				num_msg += (glob_comm_sz - my_comm_sz)%my_comm_sz;
			
			
			start_rank = (me_glob - me_comp == 0) ? me_glob + my_comm_sz : me_comp;
			
				
			for(i = 0; i < num_msg; i++){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Notify %d that my comp finished", start_rank+i);
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, COMP_FINALIZED_TAG, 1);
				err = MPI_Isend(NULL, 0, MPI_INT, start_rank + i, COMP_FINALIZED_TAG, gl_proc.comps[comp].comm, msg->reqs);
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_proc.comps[comp].out_msg_q);
			}
			 
		}
		 
	} else { // !gl_conf.single_mpirun_mode
		for(comp = 0; comp < gl_proc.ncomps; comp++){
			if(comp == gl_proc.my_comp)
				continue;
			int ncomm_my, ncomm_cpl; 
			MPI_Comm_size(gl_proc.comps[gl_proc.my_comp].comm, &ncomm_my);
			MPI_Comm_remote_size(gl_proc.comps[comp].comm, &ncomm_cpl);
			
			if(gl_proc.myrank < ncomm_cpl){
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, COMP_FINALIZED_TAG, 1);
				err = MPI_Isend(NULL, 0, MPI_INT, gl_proc.myrank, COMP_FINALIZED_TAG, gl_proc.comps[comp].comm, msg->reqs);
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_proc.comps[comp].out_msg_q);
			}
			
			if( (ncomm_cpl > ncomm_my) && (gl_proc.myrank == 0)){
				int i;
				for(i = gl_proc.myrank; i < ncomm_cpl; i++){
					dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, COMP_FINALIZED_TAG, 1);
					err = MPI_Isend(NULL, 0, MPI_INT, gl_proc.myrank, COMP_FINALIZED_TAG, gl_proc.comps[comp].comm, msg->reqs);
					CHECK_MPI(err);
					ENQUEUE_ITEM(msg, gl_proc.comps[comp].out_msg_q);
				}
			} 
		}
	}
    
    for(comp = 0; comp < gl_proc.ncomps; comp++){
		if(gl_proc.comps[comp].out_msg_q == NULL)
			continue;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Finalize message queue for comp %s", gl_proc.comps[comp].name);
		while(gl_proc.comps[comp].out_msg_q != NULL)
			progress_comm(0);	
	}
	
	gl_proc.comps[gl_proc.my_comp].finalized = 1;
    
    finalize_comp_comm();
    print_stats();
    //destroy inrracomp communicator
    err = MPI_Comm_free(&gl_proc.comps[gl_proc.my_comp].comm);
    CHECK_MPI(err);
    
	for(comp = 0; comp < gl_proc.ncomps; comp++)
		assert(gl_proc.comps[comp].out_msg_q == NULL);

    clean_config();

    if(gl_proc.msgbuf != NULL)
        dtf_free(gl_proc.msgbuf, gl_proc.conf.data_msg_size_limit);
    
	DTF_DBG(VERBOSE_DBG_LEVEL,"DTF: finalized");
    lib_initialized = 0;
    fflush(stdout);
    fflush(stderr);
    return 0;
}


/*
 * Process will block until all active data transfers are completed.
 * */
 //TODO is this still working?
_EXTERN_C_ int dtf_transfer_all_files()
{
	double t_start;
	if(!lib_initialized) return 0;
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Start transfer_complete_all");
	t_start = MPI_Wtime();
	gl_proc.stats_info.t_idle = MPI_Wtime();
	match_ioreqs_all_files();
	gl_proc.stats_info.transfer_time += MPI_Wtime() - t_start;
	gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
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
  
    fbuf = find_file_buffer(gl_proc.filebuf_list, filename, ncid);
    if(fbuf == NULL){

        if( (filename != NULL) && (strlen(filename) == 0) )
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file (%s) with ncid %d is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        else
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Matching ignored.", filename, ncid);
        return 0;
    }
	 
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    
	if(!fbuf->is_transferring) return 0;
    
    gl_proc.stats_info.t_idle = MPI_Wtime();
    
    match_ioreqs(fbuf);
     
	fname_pattern_t *pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions)
		delete_file_buffer(fbuf);
    
    gl_proc.stats_info.transfer_time += MPI_Wtime() - t_start;
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Exit transfer");
    return 0;
}

/*Supposed to be called by the writer process.
  Used to match against several dtf_transfer functions on the reader side*/
_EXTERN_C_ void dtf_transfer_multiple(const char* filename, int ncid)
{
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s ncid %d is not treated by DTF( \
                not in configuration file). Explicit matching ignored.", filename, ncid);
        return;
    }
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(!fbuf->is_transferring) return;
    
	double t_start = MPI_Wtime();
    if( gl_proc.my_comp != fbuf->writer_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_transfer_multiple can only be called by writer component.");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching multiple");
   
   	gl_proc.stats_info.t_idle = MPI_Wtime();
   
    fbuf->done_multiple_flag = 0;
    while(!fbuf->done_multiple_flag){
        match_ioreqs(fbuf);
        
        /*If the other component has started finalizing, then just complete this transfer*/
        if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
			((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
				fbuf->done_multiple_flag = 1;
		}
	}
    //reset
    fbuf->done_multiple_flag = 0;
    
	fname_pattern_t *pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions)
		delete_file_buffer(fbuf);
    
    gl_proc.stats_info.transfer_time += MPI_Wtime() - t_start;
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finish matching multiple");
    return;
}

//TODO there was a bug with timing in receiving match complete notifications in scale-letkf
/*Used by reader to notify writer that it can complete dtf_match_multiple*/
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid)
{
	fname_pattern_t *pat;
    double t_start = MPI_Wtime();
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). matching ignored.", filename, ncid);
        return;
    }

    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    
    if(fbuf->is_write_only) return;
    
    if(fbuf->reader_id != gl_proc.my_comp){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_complete_multiple for file %s can only be called by reader", filename);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    
    MPI_Barrier(fbuf->comm);
    
    gl_proc.stats_info.t_idle = MPI_Wtime();
    
    if( ((fbuf->writer_id == gl_proc.my_comp) && !gl_proc.comps[fbuf->reader_id].finalized)  ||
		((fbuf->reader_id == gl_proc.my_comp) && !gl_proc.comps[fbuf->writer_id].finalized) )			
		notify_complete_multiple(fbuf);
    
    pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(fbuf->session_cnt == pat->num_sessions && !fbuf->is_transferring)
		delete_file_buffer(fbuf);
    
    gl_proc.stats_info.transfer_time += MPI_Wtime() - t_start;
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
}

/*User controled timers*/
_EXTERN_C_ void dtf_time_start()
{
    if(!lib_initialized) return;
  
    gl_proc.stats_info.user_timer_start = MPI_Wtime();
    DTF_DBG(VERBOSE_ERROR_LEVEL, "user_time start");
}

_EXTERN_C_ void dtf_time_end()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_proc.stats_info.user_timer_start;
    gl_proc.stats_info.user_timer_accum += tt;
    gl_proc.stats_info.user_timer_start = 0;
    
	DTF_DBG(VERBOSE_ERROR_LEVEL, "user_time end  %.6f", tt);
    
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
