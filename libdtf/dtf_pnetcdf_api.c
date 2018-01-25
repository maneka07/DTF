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
#include "dtf_init_finalize.h"
#include "dtf_util.h"
#include "dtf_nbuf_io.h"
#include "dtf_req_match.h"


extern int lib_initialized;

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

    if(!lib_initialized) return;

    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf != NULL){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error creating file: a file with the same name (%s) has been created/opened before and has not been closed yet. Abort.", filename);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}

    fname_pattern_t *pat = gl_fname_ptrns;
	while(pat != NULL){
		if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls) && !pat->ignore_io){
			fbuf = create_file_buffer(pat, filename);
			/*Because this component creates the file, we assume that it's the writer*/
			fbuf->omode = DTF_WRITE; 
			fbuf->writer_id = gl_my_comp_id;
			fbuf->reader_id = (gl_my_comp_id == pat->comp1) ? pat->comp2 : pat->comp1;
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
	init_req_match_masters(comm, fbuf->my_mst_info);
	fbuf->root_writer = fbuf->my_mst_info->masters[0];
	
	 /*In scale-letkf we assume completely mirrorer file handling. 
	 * Hence, will simply duplicate the master structre*/
	 if(gl_scale){
		assert(fbuf->cpl_mst_info->nmasters == 0);
		fbuf->cpl_mst_info->nmasters = fbuf->my_mst_info->nmasters;
		fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
		assert(fbuf->cpl_mst_info->masters != NULL);
		memcpy(fbuf->cpl_mst_info->masters, fbuf->my_mst_info->masters, fbuf->cpl_mst_info->nmasters*sizeof(int));
		fbuf->root_reader = fbuf->cpl_mst_info->masters[0];
	 }
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Root master for file %s (ncid %d) is %d", filename, ncid, fbuf->root_writer);

    if(gl_my_rank == fbuf->root_writer && !gl_scale){
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
	
	if(fbuf->iomode == DTF_IO_MODE_FILE){
        if(fbuf->root_writer == gl_my_rank)
			fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;

        //scale-letkf
		if(strstr(fbuf->file_path, "hist.d")!=NULL)
			gl_stats.st_mtch_hist = MPI_Wtime()-gl_stats.walltime;
		else if(strstr(fbuf->file_path, "anal.d")!=NULL)
			gl_stats.st_mtch_rest = MPI_Wtime()-gl_stats.walltime;
		//~ char *s = getenv("DTF_SCALE");
		//~ if(s != NULL)
			//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp open file %s", fbuf->file_path);
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Exit create");
}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_open(const char *filename, int omode, MPI_Comm comm)
{
    int nranks;
    file_buffer_t *fbuf;
	
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s", filename);

    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);

    if(fbuf == NULL){
		fname_pattern_t *pat = gl_fname_ptrns;
		while(pat != NULL){
			if(match_ptrn(pat->fname, filename, pat->excl_fnames, pat->nexcls) && !pat->ignore_io){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
				fbuf = create_file_buffer(pat, filename);
				/*Assume that this component is the reader*/
				fbuf->reader_id = gl_my_comp_id;
				fbuf->writer_id = (gl_my_comp_id == pat->comp1) ? pat->comp2 : pat->comp1;
				assert( (omode & NC_WRITE) == 0);
				fbuf->omode = DTF_READ;
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

    if(fbuf->my_mst_info->nmasters == 0){ /*Open for the first time*/
        init_req_match_masters(comm, fbuf->my_mst_info);
        
        /*In scale-letkf we assume completely mirrorer file handling. 
         * Hence, will simply duplicate the master structre*/
         if(gl_scale){
			assert(fbuf->cpl_mst_info->nmasters == 0);
			fbuf->cpl_mst_info->nmasters = fbuf->my_mst_info->nmasters;
			fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
			assert(fbuf->cpl_mst_info->masters != NULL);
			memcpy(fbuf->cpl_mst_info->masters, fbuf->my_mst_info->masters, fbuf->cpl_mst_info->nmasters*sizeof(int));
			fbuf->root_writer = fbuf->cpl_mst_info->masters[0];
		 }
    } else {
		/*do a simple check that the same set of processes as before opened the file*/
		int rank;
		MPI_Comm_rank(comm, &rank);
		if(rank == 0)
			assert(gl_my_rank == fbuf->my_mst_info->masters[0]);			
	}

	/*Check if the file was opened before and if we need to 
	 * complete the file-ready notification*/
	 if(fbuf->fready_notify_flag == RDR_NOTIF_POSTED){
		assert(fbuf->root_writer == gl_my_rank);
		while(fbuf->fready_notify_flag != DTF_UNDEFINED)
			progress_io_matching();
	 }

	if( omode & NC_WRITE )
		fbuf->omode = DTF_WRITE;
	else // NC_NOWRITE 
		fbuf->omode = DTF_READ;

	/*Set who's the reader and writer component in this session*/
	int cpl_cmp = (fbuf->reader_id == gl_my_comp_id) ? fbuf->writer_id : fbuf->reader_id; 
	int cpl_root = (fbuf->reader_id == gl_my_comp_id) ? fbuf->root_writer : fbuf->root_reader;
	
	if(fbuf->omode == DTF_READ){

	    /*If the component was the writer in the previous session, the root
		 * process must have the master info of the coupled component. It
		 * will broadcast the info.*/
		if(fbuf->writer_id == gl_my_comp_id){
			int rank, err;
			
			if(!gl_scale){
				MPI_Comm_rank(comm, &rank);
				DTF_DBG(VERBOSE_DBG_LEVEL, "Broadcast info about the other component");
				
				//root broadcasts master info to others
				if(rank == 0){
					assert(fbuf->cpl_mst_info->nmasters>0);
				}
				err = MPI_Bcast(&(fbuf->cpl_mst_info->nmasters), 1, MPI_INT, 0, comm);
				CHECK_MPI(err);
				if(rank != 0){
					assert(fbuf->cpl_mst_info->masters== NULL);
					fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
					assert(fbuf->cpl_mst_info->masters != NULL);	
				}
				
				err = MPI_Bcast(fbuf->cpl_mst_info->masters, fbuf->cpl_mst_info->nmasters, MPI_INT, 0, comm);
				CHECK_MPI(err);
			}
			
			fbuf->root_writer = fbuf->cpl_mst_info->masters[0];
			
			//~ for(i = 0; i < fbuf->cpl_mst_info->nmasters; i++)
				//~ DTF_DBG(VERBOSE_DBG_LEVEL, "mst list %d",fbuf->cpl_mst_info->masters[i]); 
		} else 
			fbuf->root_writer = cpl_root;
		 
		fbuf->reader_id = gl_my_comp_id;
		fbuf->root_reader = fbuf->my_mst_info->masters[0];
		fbuf->writer_id = cpl_cmp;
	} else { /*fbuf->omode == DTF_WRITE*/			
		fbuf->writer_id = gl_my_comp_id;
		fbuf->root_writer = fbuf->my_mst_info->masters[0];
		fbuf->reader_id = cpl_cmp;	
		fbuf->root_reader = cpl_root;
	}
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Writer %d (root %d), reader %d (root %d), omode %d", fbuf->writer_id, fbuf->root_writer, fbuf->reader_id, fbuf->root_reader, omode);
    
    open_file(fbuf, comm);
}

_EXTERN_C_ void dtf_enddef(const char *filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;

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

    close_file(fbuf);
}

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
    if(rw_flag==DTF_WRITE && fbuf->reader_id == gl_my_comp_id){
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
	int i;
	if(!lib_initialized) return;
	if(!gl_conf.log_ioreqs) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_FILE) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "log ioreq for var %d", varid);
    for(i=0; i < ndims; i++)
		DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", start[i], count[i]);
    log_ioreq(fbuf, varid, ndims, start, count, dtype, buf, rw_flag);

}



/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/

_EXTERN_C_ int dtf_io_mode(const char* filename)
{
    if(!lib_initialized) return DTF_UNDEFINED;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        return DTF_UNDEFINED;
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
   // DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat %.6f", tt);
 
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
