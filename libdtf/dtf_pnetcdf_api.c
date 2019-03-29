/*
* Interfaces to be used inside the PnetCDF library
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
#include "dtf_io_pattern.h"


extern int lib_initialized;

_EXTERN_C_ void dtf_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    
    double t_start = MPI_Wtime();
    
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;

    if(hdr_sz == 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Header size for file %s is zero", filename);
        return;
    }

    write_hdr(fbuf, hdr_sz, header);
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
}

_EXTERN_C_ MPI_Offset dtf_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
	MPI_Offset ret;
    if(!lib_initialized) return 0;

    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    double t_start = MPI_Wtime();
    
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
	
    ret = read_hdr_chunk(fbuf, offset, chunk_sz, chunk);
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
    return ret;
}

_EXTERN_C_ void dtf_create(const char *filename, MPI_Comm comm)
{
    file_buffer_t *fbuf;
	double t_start = MPI_Wtime();
	
    if(!lib_initialized) return;
    
	gl_proc.stats_info.t_idle = MPI_Wtime();

    fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf != NULL){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error creating file: a file with the same name (%s) has been created/opened before and has not been closed yet. Abort.", filename);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	
    fname_pattern_t *pat = find_fname_pattern(filename);
    if(pat != NULL){
		fbuf = create_file_buffer(pat, filename, comm);
		/*Because this component creates the file, we assume that it's the writer*/
		//fbuf->omode = DTF_WRITE; 
		fbuf->writer_id = gl_proc.my_comp;
		fbuf->reader_id = (gl_proc.my_comp == pat->comp1) ? pat->comp2 : pat->comp1;
	}

	if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Created file %s. File is not treated by DTF", filename);
        return;
    } else {
        DTF_DBG(VERBOSE_DBG_LEVEL, "Creating file %s", filename);
    }
	
    fbuf->comm = comm;
    
	fbuf->root_writer = fbuf->my_mst_info->masters[0];
	if(pat->replay_io && (pat->wrt_recorded == IO_PATTERN_RECORDING)){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: Creating file %s but have not finished recording the I/O pattern for \
				a file that matches the same file pattern. Did you forget to close the old file?", filename);
		pat->wrt_recorded = IO_PATTERN_RECORDED;
	} else if (pat->replay_io && (pat->wrt_recorded == DTF_UNDEFINED))
		pat->wrt_recorded = IO_PATTERN_RECORDING;
			 
	 if(pat->replay_io && (pat->finfo_sz > 0)){
		int offt = MAX_FILE_NAME;
		void *rbuf = pat->finfo;
		//parse mst info of the other component
		assert(fbuf->cpl_mst_info->nmasters == 0);
		memcpy(&(fbuf->cpl_mst_info->comm_sz), (unsigned char*)rbuf+offt, sizeof(int));
		offt+=sizeof(int);
		memcpy(&(fbuf->cpl_mst_info->nmasters), (unsigned char*)rbuf+offt, sizeof(int));
		offt+=sizeof(int);
		assert(fbuf->cpl_mst_info->nmasters > 0);
		fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
		memcpy(fbuf->cpl_mst_info->masters, (unsigned char*)rbuf+offt, fbuf->cpl_mst_info->nmasters*sizeof(int));
		fbuf->root_reader = fbuf->cpl_mst_info->masters[0];
		fbuf->cpl_info_shared = 1;		
	 }
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Root master for file %s is %d", filename, fbuf->root_writer);

    if(!pat->mirror_io_root && (gl_proc.myrank == fbuf->root_writer) && (fbuf->cpl_mst_info->nmasters == 0)){
		FILE *rootf;
		char outfname[MAX_FILE_NAME+5];
		strcpy(outfname, filename);
		strcat(outfname, ".root");
		DTF_DBG(VERBOSE_DBG_LEVEL, "Will write root to file %s", outfname);
		//we assume that there are no race conditions
		rootf = fopen(outfname, "wb");
		assert(rootf != NULL);
		if(!fwrite(&fbuf->root_writer, sizeof(int), 1, rootf))
			assert(0);
		fclose(rootf);
	}
	
	if(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->root_writer == gl_proc.myrank)
			fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;
    
    
	DTF_DBG(VERBOSE_DBG_LEVEL, "Writer %s (root %d), reader %s (root %d)", gl_proc.comps[fbuf->writer_id].name, fbuf->root_writer, gl_proc.comps[fbuf->reader_id].name, fbuf->root_reader);	
    DTF_DBG(VERBOSE_DBG_LEVEL, "Exit create");
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_open(const char *filename, int omode, MPI_Comm comm)
{
    int nranks, cpl_cmp, cpl_root;
    file_buffer_t *fbuf;
	fname_pattern_t *pat = NULL;
	
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s", filename);
    
	gl_proc.stats_info.t_idle = MPI_Wtime();

    fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);

    if(fbuf == NULL){
		pat = find_fname_pattern(filename);
		if(pat != NULL){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
			fbuf = create_file_buffer(pat, filename, comm);	
		}
	} else {
		assert(fbuf->comm == comm); //opened in the same communicator as before
	}

    if(fbuf == NULL) {
        DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s. File is not treated by DTF", filename);
        return;
    }
    
    fbuf->is_transferring = 0;  //reset

	progress_comm(1);

    MPI_Comm_size(comm, &nranks);

    if(comm == gl_proc.comps[gl_proc.my_comp].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in component's communicator (%d nprocs)", nranks);
    else
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in subcommunicator (%d nprocs)", nranks);

    

	//~ if(fbuf->my_mst_info->nmasters == 0){ 
		/*Open for the first time*/        
		/*In scale-letkf we assume completely mirrorer file handling. 
		 * Hence, will simply duplicate the master structre*/
	
    //~ } else {
		//~ /*do a simple check that the same set of processes as before opened the file*/
		//~ int rank;
		//~ MPI_Comm_rank(comm, &rank);
		//~ if(rank == 0)
			//~ assert(gl_proc.myrank == fbuf->my_mst_info->masters[0]);			
	//~ }

	/*Check if the file was opened before and if we need to 
	 * complete the file-ready notification*/
	 if(fbuf->fready_notify_flag == RDR_NOTIF_POSTED){
		assert(fbuf->root_writer == gl_proc.myrank);
		while(fbuf->fready_notify_flag != RDR_NOTIFIED)
			progress_comm(0);
		//reset
		fbuf->fready_notify_flag = DTF_UNDEFINED;
	 }	   

	if(pat == NULL){
		pat = find_fname_pattern(filename);
		assert(pat != NULL);
		DTF_DBG(VERBOSE_DBG_LEVEL, "Pat comp1 %d, comp2 %d, fname %s, replay %d", pat->comp1, pat->comp2, pat->fname, pat->replay_io);
	}
	 
	if(pat->replay_io){
		 if(pat->wrt_recorded == IO_PATTERN_RECORDING || pat->rdr_recorded == IO_PATTERN_RECORDING)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: Opening file %s but have not finished recording the I/O pattern for \
			a file that matches the same file pattern. Did you forget to close the old file?", filename);
		 if(pat->wrt_recorded == IO_PATTERN_RECORDING)
			pat->wrt_recorded = IO_PATTERN_RECORDED;
		 if(pat->rdr_recorded == IO_PATTERN_RECORDING)
			pat->rdr_recorded = IO_PATTERN_RECORDED;
		
		//Get the file info (header, list of vars etc. and info about the other component) 
		//saved from previous session. Assume that this information is unchanged in this session.
		if(pat->finfo_sz > 0 && fbuf->header == NULL)
			unpack_file_info(pat->finfo_sz, pat->finfo, fbuf);
	}
	
	/*Set who's the reader and writer component in this session*/
	cpl_cmp = (pat->comp1 == gl_proc.my_comp) ? pat->comp2 : pat->comp1; 
	
	/*if this is not the first transfer session, the root from the other coponent is 
	 * known from previous session*/
	if(fbuf->cpl_mst_info->nmasters > 0)
		cpl_root = fbuf->cpl_mst_info->masters[0];
	else 
		cpl_root = (fbuf->reader_id == gl_proc.my_comp) ? fbuf->root_writer : fbuf->root_reader;
	
	if( omode & NC_WRITE && !fbuf->is_write_only){
		//fbuf->omode = DTF_WRITE;
		fbuf->writer_id = gl_proc.my_comp;
		fbuf->reader_id = cpl_cmp;
		fbuf->root_writer = fbuf->my_mst_info->masters[0];
		fbuf->root_reader = cpl_root;
		
		if(pat->replay_io && (pat->wrt_recorded == DTF_UNDEFINED))
			pat->wrt_recorded = IO_PATTERN_RECORDING;
			
	} else { // NC_NOWRITE  
		
		//fbuf->omode = DTF_READ;
		fbuf->writer_id = cpl_cmp;
		fbuf->reader_id = gl_proc.my_comp;
		fbuf->root_writer = cpl_root;
		fbuf->root_reader = fbuf->my_mst_info->masters[0];
	
		if(pat->replay_io && (pat->rdr_recorded == DTF_UNDEFINED))
			pat->rdr_recorded = IO_PATTERN_RECORDING;
		
		if(pat->write_only)
			DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is a write only file, consider as if opened for reading", filename);
	}
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Writer %s (root %d), reader %s (root %d), omode %d", gl_proc.comps[fbuf->writer_id].name, fbuf->root_writer, gl_proc.comps[fbuf->reader_id].name, fbuf->root_reader, omode);	  

	open_file(fbuf, comm);
    //reset 
    gl_proc.stats_info.t_idle = MPI_Wtime();
    DTF_DBG(VERBOSE_DBG_LEVEL,"Exit open");
}

_EXTERN_C_ void dtf_enddef(const char *filename)
{
    if(!lib_initialized) return;
    double t_start = MPI_Wtime();
    
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    fbuf->is_defined = 1;
    progress_comm(1);
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
}

_EXTERN_C_ void dtf_set_ncid(const char *filename, int ncid)
{
    if(!lib_initialized) return;
	double t_start = MPI_Wtime();
	
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is not treated by DTF. Will not set ncid", filename);
        return;
    }
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Set ncid of file %s to %d (previos value %d)", filename, ncid, fbuf->ncid);
    fbuf->ncid = ncid;
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
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
	fname_pattern_t *pat;
	
    if(!lib_initialized) return;
    double t_start = MPI_Wtime();
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
   
    if(fbuf == NULL) return;
    
	gl_proc.stats_info.t_idle = MPI_Wtime();
		
	progress_comm(1);

	if(fbuf->ioreq_log != NULL){ // && strstr(fbuf->file_path, "anal")!=NULL){
		double check;
		int i;
		
		io_req_log_t *req = fbuf->ioreq_log;
		while(req != NULL){
			
			DTF_DBG(VERBOSE_DBG_LEVEL, "Req for var %d", req->var_id);
			
			for(i=0; i <req->ndims; i++)
				DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", req->start[i], req->count[i]);
			
			check = compute_checksum(req->user_buf, req->ndims, req->count, req->dtype);
			DTF_DBG(VERBOSE_DBG_LEVEL, "Checksum %.3f", check);
			req = req->next;
		}
	}

	close_file(fbuf);
	
	fbuf->cur_transfer_epoch = 0; //reset
	fbuf->session_cnt++;
	
	pat = find_fname_pattern(filename);			
	assert(pat != NULL);
	if(pat->replay_io){
		//need to confirm that io pattern for this eposh has been recorded
		//if file is closed before transfer function is called then 
		//pattern is not recorded yet
		if(pat->rdr_recorded == IO_PATTERN_RECORDING)
			pat->rdr_recorded = IO_PATTERN_RECORDED;
		else if(pat->wrt_recorded == IO_PATTERN_RECORDING) 
			pat->wrt_recorded = IO_PATTERN_RECORDED;
	}
	
	if(fbuf->session_cnt == pat->num_sessions){
		if( (fbuf->iomode == DTF_IO_MODE_MEMORY && !fbuf->is_transferring) || 
			(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->root_writer == gl_proc.myrank && fbuf->fready_notify_flag == RDR_NOTIFIED) ||
			(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->root_writer != gl_proc.myrank) || 
			(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->reader_id == gl_proc.my_comp)){
				if(gl_proc.stats_info.nwreqs > 0 || gl_proc.stats_info.nrreqs > 0 )
					DTF_DBG(VERBOSE_DBG_LEVEL, "%lu rreqs and %lu wreqs, %lu my own reqs", gl_proc.stats_info.nrreqs	,  gl_proc.stats_info.nwreqs, gl_proc.stats_info.nioreqs);
				delete_file_buffer(fbuf);
			}
	
	}
	
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
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
    
    double t_start = MPI_Wtime();
    
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
  
    
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
	
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
    
	if(fbuf->is_write_only){
		DTF_DBG(VERBOSE_DBG_LEVEL, "ignore_io or write_only options are set for this file. I/O for this file will be ignored");
		int el_sz, i;
        dtf_var_t *var = fbuf->vars[varid];
        MPI_Type_size(dtype, &el_sz);
		ret = 1;
		for(i = 0; i < var->ndims; i++)
			ret *= count[i];
		ret *= el_sz;
		return ret;
	}
    
    if(rw_flag == DTF_WRITE && fbuf->reader_id == gl_proc.my_comp){
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
    
    if(rw_flag == DTF_READ){
        if(fbuf->reader_id == gl_proc.my_comp){
          assert(fbuf->is_ready);
        } else{
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: writer process tries to read file %s (var %d)", fbuf->file_path, varid);
            //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(imap != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing mapped vars is not impelemented yet. Ignore.");
        assert(0);
    }

    if(stride != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing vars at a stride is not impelemented yet. Ignore.");
       assert(0);
    }

    ret = read_write_var(fbuf, varid, start, count, dtype, buf, rw_flag);
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
    fbuf->is_transferring = 1;
    return ret;
}


_EXTERN_C_ void dtf_log_ioreq(const char *filename,
                                          int varid, int ndims,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag)
{
	int i, type_sz;
	if(!lib_initialized) return;
	if(!gl_proc.conf.log_ioreqs) return;
	double t_start = MPI_Wtime();
	
    file_buffer_t *fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    //if(fbuf->iomode != DTF_IO_MODE_FILE) return;
    MPI_Type_size(dtype, &type_sz);
    DTF_DBG(VERBOSE_DBG_LEVEL, "log ioreq for var %d, type sz %d, rw %d", varid, type_sz, rw_flag);
    for(i=0; i < ndims; i++)
		DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", start[i], count[i]);
    log_ioreq(fbuf, varid, ndims, start, count, dtype, buf, rw_flag);
	gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
}



/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/

_EXTERN_C_ int dtf_io_mode(const char* filename)
{
    if(!lib_initialized) return DTF_UNDEFINED;
    
    fname_pattern_t * pat = find_fname_pattern(filename);	
    if(pat == NULL){
        return DTF_UNDEFINED;
    }
    return pat->iomode;
}

_EXTERN_C_ int dtf_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i, ret;
    double t_start = MPI_Wtime();
    
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
   
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
    gl_proc.stats_info.dtf_time += MPI_Wtime() - t_start;
    return ret;
}

/*Library controled timers*/
_EXTERN_C_ void dtf_tstart()
{
    if(!lib_initialized) return;
    
    if(gl_proc.stats_info.timer_start != 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: user timer was started at %.3f and not finished.",
                gl_proc.stats_info.timer_start - gl_proc.walltime);

    gl_proc.stats_info.timer_start = MPI_Wtime();
}
_EXTERN_C_ void dtf_tend()
{
    double tt;
    if(!lib_initialized) return;
    tt = MPI_Wtime() - gl_proc.stats_info.timer_start;
    gl_proc.stats_info.timer_accum += tt;
    gl_proc.stats_info.timer_start = 0;
   // DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat %.6f", tt);
 
}

_EXTERN_C_ MPI_File *dtf_get_tmpfile()
{
	return &gl_proc.tmpfile;
}
