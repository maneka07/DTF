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
	return 0;
    //~ int ret;
    //~ if(!lib_initialized) return 0;
    //~ file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    //~ if(fbuf == NULL) return 0;
    //~ if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    //~ if(fbuf->ignore_io) return 0;

     //~ //Never checked if this function works correctly
     //~ assert(0);

    //~ if(fbuf->is_matching_flag){
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_ioreqs is called for file %s, but a matching process has already started before.", fbuf->file_path);
        //~ MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    //~ }
    //~ fbuf->is_matching_flag = 1;
    //~ ret = match_ioreqs(fbuf, 0);
    //~ fbuf->is_matching_flag = 0;

    //~ return ret;
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
    DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat %.6f", tt);
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
