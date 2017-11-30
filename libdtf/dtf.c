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

static int lib_initialized=0;

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
struct dtf_config gl_conf;                 /*Framework settings*/
struct stats gl_stats;
char *gl_my_comp_name = NULL;
void* gl_msg_buf = NULL;

/*
************************************User API**********************************************
*/
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
    DTF_DBG(VERBOSE_DBG_LEVEL,"DTF: finalize");
	
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
    //TODO send notification to root 0 when file is closed so that it's deleted from the list
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
            DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file (%s) with ncid %d is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
        else
            DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
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
    if(fbuf->ignore_io) return 0;
    if(!fbuf->explicit_match){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: calling dtf_match_io but explicit match for file %s not enabled. Ignored.", filename);
        return 0;
    }

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
	file_buffer_t *fbuf;
	if(!lib_initialized) return;
	
	if(filename == NULL || strlen(filename)==0){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: undefined file name in dtf_skip_match call");
		return;
	}
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Call skip match for %s", filename);
	
	fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
		//If fbuf doesn't exist this file hasn't been opened/created. 
		//However, we still need to check if I/O matching is enabled for 
		//this file and, if it is, notify the other component that we skip.
		fname_pattern_t *pat = gl_fname_ptrns;
		while(pat != NULL){
			if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls)){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
				break;
			}
			pat = pat->next;
		} 
		
		if(pat != NULL){
			if(gl_my_comp_id == pat->wrt){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_skip_match for file %s should not be called by the writer component. Ignore. ", filename);
				return;
			}
			assert(pat->rdr == gl_my_comp_id);
			
			if(pat->iomode == DTF_IO_MODE_MEMORY){	
				assert(comm != MPI_COMM_NULL);
				skip_match(NULL, filename, comm, pat->wrt);
			} 	
		}
       
    } else if(fbuf->iomode == DTF_IO_MODE_MEMORY)
		skip_match(fbuf, filename, comm, -1);
	
}

/*Supposed to be called by the writer process.
  Used to match against several dtf_match_io functions on the reader side*/
_EXTERN_C_ void dtf_match_multiple(int ncid)
{
	//TODO figure out with cleaning iodb and notifications
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: ncid %d is not treated by DTF( \
                not in configuration file). Explicit matching ignored.", ncid);
        return;
    }
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(fbuf->ignore_io) return;
    if(!fbuf->explicit_match){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: calling dtf_match_multiple but explicit match for file \
                %s not enabled. Ignored.", fbuf->file_path);
        return;
    }

    if( gl_my_comp_id != fbuf->writer_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_match_multiple can only be called by writer component. Ignoring.");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching multiple");
    dtf_tstart();
    fbuf->done_match_multiple_flag = 0;
    while(!fbuf->done_match_multiple_flag)
        match_ioreqs(fbuf, 0);
    //reset
    fbuf->done_match_multiple_flag = 0;
    dtf_tend();
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finish matching multiple");
    return;
}

//TODO remove this?
/*Used by reader to notify writer that it can complete dtf_match_multiple*/
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid)
{
    double t_start;
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
        return;
    }

    t_start = MPI_Wtime();
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(fbuf->reader_id != gl_my_comp_id){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: dtf_complete_multiple can only be called by reader");
        return;
    }
    dtf_tstart();
    //MPI_Barrier(fbuf->comm);
    notify_complete_multiple(fbuf);
    dtf_tend();

    gl_stats.accum_match_time += MPI_Wtime() - t_start;
}

/*
*******************************Interfaces to be used by PnetCDF***************************************************
*/
_EXTERN_C_ void dtf_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;


    if(hdr_sz == 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Header size for file %s is zero", filename);
        return;
    }
    double t_st = MPI_Wtime();
    write_hdr(fbuf, hdr_sz, header);
    gl_stats.accum_hdr_time += MPI_Wtime() - t_st;

}

_EXTERN_C_ MPI_Offset dtf_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
	MPI_Offset ret;
    if(!lib_initialized) return 0;
 
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    double t_st = MPI_Wtime();
    ret = read_hdr_chunk(fbuf, offset, chunk_sz, chunk);
    gl_stats.accum_hdr_time += MPI_Wtime() - t_st;
    return ret;
}

_EXTERN_C_ void dtf_create(const char *filename, MPI_Comm comm, int ncid)
{   
    file_buffer_t *fbuf;
    fname_pattern_t *pat;
    
    if(!lib_initialized) return;
  
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf != NULL){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error creating file: a file with the same name (%s) has been created/opened before and has not been closed yet. Abort.", filename);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
   
    pat = gl_fname_ptrns;
	while(pat != NULL){
		if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls)){
			fbuf = create_file_buffer(pat, filename);
			break;
		}
		pat = pat->next;
	} 
    
	if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Creating file %s. File is not treated by DTF", filename);
        return;
    } else {
        DTF_DBG(VERBOSE_DBG_LEVEL, "Created file %s (ncid %d)", filename, ncid);
    }
    
    fbuf->ncid = ncid;
    fbuf->comm = comm;
    int root_mst = gl_my_rank;
    int err = MPI_Bcast(&root_mst, 1, MPI_INT, 0, comm);
    CHECK_MPI(err);
    fbuf->root_writer = root_mst;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Root master for file %s (ncid %d) is %d", filename, ncid, fbuf->root_writer);
    //~ if(gl_my_rank == fbuf->root_writer){
        //~ /*Notify global rank 0 that I am the root master for this file*/
        //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Will notify global rank 0 that I am master");
        //~ dtf_msg_t *msg = new_dtf_msg(NULL, 0, ROOT_MST_TAG);
        //~ err = MPI_Isend(fbuf->file_path, (int)(MAX_FILE_NAME), MPI_CHAR, 0, ROOT_MST_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
        //~ CHECK_MPI(err);
        //~ ENQUEUE_ITEM(msg, gl_msg_q);
    //~ }
    
    if(gl_my_rank == fbuf->root_writer){
		int i;
		FILE *rootf;
		char *glob = getenv("DTF_GLOBAL_PATH");
		assert(glob != NULL);
		char outfname[MAX_FILE_NAME*2];
		char fname[MAX_FILE_NAME];
		strcpy(fname, filename);
		for(i = 0; i <strlen(filename); i++)
			if(fname[i]=='/' || fname[i]=='\\')
				fname[i]='_';
		strcpy(outfname, glob);
		strcat(outfname, "/");
		strcat(outfname, fname);
		strcat(outfname, ".root");
		DTF_DBG(VERBOSE_DBG_LEVEL, "Will write root to file %s", outfname);
		//we assume that there are no race conditions
		rootf = fopen(outfname, "wb");
		assert(rootf != NULL);
		if(!fwrite(&fbuf->root_writer, sizeof(int), 1, rootf))
			assert(0);
		fclose(rootf);
		
	}
    

    DTF_DBG(VERBOSE_DBG_LEVEL, "Init masters");
    init_req_match_masters(comm, fbuf->mst_info);

    if(fbuf->iomode == DTF_IO_MODE_MEMORY){
        if(fbuf->mst_info->is_master_flag){
            int nranks;
            MPI_Comm_size(comm, &nranks);

            assert(fbuf->mst_info->iodb == NULL);
            init_iodb(fbuf);
            fbuf->mst_info->nranks_opened = (unsigned int)nranks;
        }

    } else if(fbuf->iomode == DTF_IO_MODE_FILE){
        fbuf->fready_notify_flag = RDR_HASNT_OPENED;

        /*Create symlink to this file (needed for SCALE-LETKF since
          there is no way to execute the script to perform all the file
          movement in the middle of the execution)*/
         if(fbuf->slink_name!=NULL && fbuf->root_writer==gl_my_rank){
            int err;
            size_t slen1, slen2;
            char *dir = NULL;
            char wdir[MAX_FILE_NAME]="\0";
            char slink[MAX_FILE_NAME]="\0";
            char origfile[MAX_FILE_NAME]="\0";

            dir = getcwd(wdir, MAX_FILE_NAME);
            assert(dir != NULL);
            slen1 = strlen(wdir);
            slen2 = strlen(fbuf->slink_name);
            if(slen1+slen2+1 > MAX_FILE_NAME){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: symlink name of length %ld exceeds max filename of %d",slen1+slen2+1, MAX_FILE_NAME);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }

            sprintf(slink, "%s/%s", wdir, fbuf->slink_name);
            sprintf(origfile, "%s/%s", wdir, fbuf->file_path);

            DTF_DBG(VERBOSE_DBG_LEVEL, "Creating symlink %s to %s (%s)", slink, origfile, fbuf->file_path);
            err = symlink(origfile, slink);
            if(err != 0){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: error creating symlink %s to %s : %s", slink, origfile,  strerror(errno));

            }
        }
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Exit create");
}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_open(const char *filename, MPI_Comm comm)
{
    int nranks;
    file_buffer_t *fbuf;

    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s", filename);
    
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    
    if(fbuf == NULL){  
		fname_pattern_t *pat = gl_fname_ptrns;
		while(pat != NULL){
			if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls)){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
				fbuf = create_file_buffer(pat, filename);
				break;
			}
			pat = pat->next;
		} 
	}
    
    if(fbuf == NULL) {
        DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s. File is not treated by DTF", filename);
        return;
    }
    

    if(fbuf->comm == MPI_COMM_NULL)
        fbuf->comm = comm;
    else
        assert(fbuf->comm == comm);
    MPI_Comm_size(comm, &nranks);
    if(comm == gl_comps[gl_my_comp_id].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in component's communicator (%d nprocs)", nranks);
    else
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in subcommunicator (%d nprocs)", nranks);
    open_file(fbuf, comm);
}

_EXTERN_C_ void dtf_enddef(const char *filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;


	if(fbuf->mst_info->is_master_flag){
		int i;
		assert(fbuf->mst_info->iodb != NULL);

		fbuf->mst_info->iodb->witems = dtf_malloc(fbuf->nvars*sizeof(write_db_item_t*));
		assert(fbuf->mst_info->iodb->witems != NULL);
		for(i = 0; i < fbuf->nvars; i++)
			fbuf->mst_info->iodb->witems[i] = NULL;
	}
}

_EXTERN_C_ void dtf_set_ncid(const char *filename, int ncid)
{
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is not treated by DTF. Will not set ncid", filename);
        return;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Set ncid of file %s to %d (previos value %d)", filename, ncid, fbuf->ncid);
    fbuf->ncid = ncid;
}

/**
  @brief	Called when the corresponding file is closed. If the file is opened in write mode
            (output file) we will notify the reader that the file is ready to be transfered. If it's opened in read mode
            (input file), we will free the buffer.
            The mode is defined from the dtf configuration file.
  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_close(const char* filename)
{
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File not treated by dtf");
        return;
    }
//    if(fbuf->reader_id == gl_my_comp_id)

    //MPI_Barrier(fbuf->comm);

    close_file(fbuf);
}

//TODO delete this I/O
/*called inside wait function in pnetcdf*/
_EXTERN_C_ int dtf_match_ioreqs(const char* filename)
{
    int ret;
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    if(fbuf->ignore_io) return 0;
    /*User will have to explicitly initiate matching*/
    if(fbuf->explicit_match) return 0;

     //Never checked if this function works correctly
     assert(0);

    if(fbuf->is_matching_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_ioreqs is called for file %s, but a matching process has already started before.", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    fbuf->is_matching_flag = 1;
    ret = match_ioreqs(fbuf, 0);
    fbuf->is_matching_flag = 0;

    return ret;
}



///* The user has to state whether the process needs to match all read or write requests.
//   Because the process of matching for a reader and writer is not the same. */
//_EXTERN_C_ void dtf_match_io_all(int rw_flag)
//{
//    file_buffer_t *fbuf = gl_filebuf_list;
//
//    if(!lib_initialized) return;
//    if(gl_conf.distr_mode != DISTR_MODE_REQ_MATCH) return;
//
//    if(rw_flag == DTF_READ){
//        DTF_DBG(VERBOSE_WARNING_LEVEL, "dtf_match_io_all() cannot be used in reader component. Ignoring.");
//        return;
//    }
//    /*First check that no matching is already happening*/
//    while(fbuf != NULL){
//        if(fbuf->is_matching_flag){
//            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_io_all is called, but a matching process has already started before.");
//            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
//        }
//        fbuf = fbuf->next;
//    }
//
//    match_ioreqs_all(rw_flag);
//
//    return;
//}


_EXTERN_C_ void dtf_print_data(int varid, int dtype, int ndims, MPI_Offset* count, void* data)
{
    return;
//
//    int i, nelems=1, max_print;
//    if(count == NULL)
//        return;
//    for(i = 0; i < ndims; i++)
//        nelems*=count[i];
//
//    if(nelems < 20)
//        max_print = nelems;
//    else
//        max_print = 20;
//    DTF_DBG(VERBOSE_ERROR_LEVEL, "Data for var %d:", varid);
//    for(i = 0; i < max_print; i++)
//        if(dtype == 0)
//            printf("%.3f\t", ((float*)data)[i]);
//        else if(dtype == 1)
//            printf("%.3f\t", ((double*)data)[i]);
//    printf("\n");
}


_EXTERN_C_ void dtf_print(const char *str)
{
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_ERROR_LEVEL, "%s", str);
}

_EXTERN_C_ MPI_Offset dtf_read_write_var(const char *filename,
                                          int varid,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          const MPI_Offset *stride,
                                          const MPI_Offset *imap,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag)
{
    MPI_Offset ret;

    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY)
        return 0;

	if(varid >= fbuf->nvars){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: variable with id %d does not exist. Abort.", varid);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}

    if(boundary_check(fbuf, varid, start, count ))
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    if( rw_flag != DTF_READ && rw_flag != DTF_WRITE){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: rw_flag value incorrect (%d)", rw_flag);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(rw_flag==DTF_WRITE && (fbuf->reader_id == gl_my_comp_id || fbuf->ignore_io)){
        int el_sz, i;
        dtf_var_t *var = fbuf->vars[varid];
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: reader component cannot write to the file %s. Ignore I/O call", fbuf->file_path);
        MPI_Type_size(dtype, &el_sz);
		ret = 1;
		for(i = 0; i < var->ndims; i++)
			ret *= count[i];
		ret *= el_sz;
		return ret;
        //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

   
        return nbuf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag);
}


_EXTERN_C_ void dtf_log_ioreq(const char *filename,
                                          int varid, int ndims,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag)
{
	if(!lib_initialized) return;
	if(!gl_conf.log_ioreqs) return;
	
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_FILE) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "log ioreq");
    log_ioreq(fbuf, varid, ndims, start, count, dtype, buf, rw_flag);
    
}

/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/

_EXTERN_C_ int dtf_io_mode(const char* filename)
{
    if(!lib_initialized) return DTF_IO_MODE_UNDEFINED;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        return DTF_IO_MODE_UNDEFINED;
    }
    return fbuf->iomode;
}

_EXTERN_C_ int dtf_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i, ret;
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;

    DTF_DBG(VERBOSE_DBG_LEVEL, "def var %d for ncid %d", varid, fbuf->ncid);

    if( (ndims > 0) && (shape[ndims - 1] == DTF_UNLIMITED))
        DTF_DBG(VERBOSE_DBG_LEVEL, "var has unlimited dimension");

    for(i = 1; i < ndims-1; i++){
        //we can support unlimited dimension if it's the first dimension
        //else abort
        if(shape[i] == DTF_UNLIMITED){
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: cannot support when unlimited dimension is not in the slowest changing dimension (dim 0). Aborting.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(dtype == MPI_DATATYPE_NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: datatype for var %d (file %s) is null. Aborting.", varid, fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    ret = def_var(fbuf, varid, ndims, dtype, shape);

    return ret;
}

/*Library controled timers*/
_EXTERN_C_ void dtf_tstart()
{
    if(!lib_initialized) return;
    //DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat start");
    if(gl_stats.timer_start != 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: user timer was started at %.3f and not finished.",
                gl_stats.timer_start - gl_stats.walltime);

    gl_stats.timer_start = MPI_Wtime();
}
_EXTERN_C_ void dtf_tend()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_stats.timer_start;
    gl_stats.timer_accum += tt;
    gl_stats.timer_start = 0;
    DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stat %.6f", tt);
 //   DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat: user time %.4f", tt);
}

/*User controled timers*/
_EXTERN_C_ void dtf_time_start()
{
    if(!lib_initialized) return;
  
    if(gl_stats.user_timer_start != 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: user timer was started at %.3f and not finished.",
                gl_stats.user_timer_start - gl_stats.walltime);

    gl_stats.user_timer_start = MPI_Wtime();
}
_EXTERN_C_ void dtf_time_end()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_stats.user_timer_start;
    gl_stats.user_timer_accum += tt;
    gl_stats.user_timer_start = 0;
    DTF_DBG(VERBOSE_ERROR_LEVEL, "user_time_stat  %.6f", tt);
 //   DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat: user time %.4f", tt);
}

/************************************************  Fortran Interfaces  *********************************************************/

_EXTERN_C_ void dtf_tstart_()
{
    dtf_tstart();
}
_EXTERN_C_ void dtf_tend_()
{
    dtf_tend();
}

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

//void dtf_match_io_all_(int *rw_flag)
//{
//    dtf_match_io_all(*rw_flag);
//}

void dtf_print_(const char *str)
{
    dtf_print(str);
}

void dtf_match_multiple_(int *ncid)
{
    dtf_match_multiple(*ncid);
}

void dtf_complete_multiple_(const char *filename, int *ncid)
{
    dtf_complete_multiple(filename, *ncid);
}

void dtf_print_data_(int *varid, int *dtype, int *ndims, MPI_Offset* count, void* data)
{
    dtf_print_data(*varid, *dtype, *ndims, count, data);
}
