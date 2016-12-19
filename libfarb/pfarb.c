#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stddef.h>

#include "pfarb_common.h"
#include "pfarb.h"
#include "pfarb_init_finalize.h"
#include "pfarb_util.h"
#include "pfarb_buf_io.h"
#include "pfarb_nbuf_io.h"
#include "pfarb_req_match.h"
//TODO !!!!!!!!!!!!!!!!!!!!!!!!implement a tree instead of list

//TODO check for dimension size and not tresspasing when writing
//TODO rename all to pfarb

int lib_initialized=0;
int gl_verbose;
int gl_my_rank;
struct farb_config gl_conf;
int frt_indexing = 0;
int matching_flag = 0;

/**
  @brief	Function to initialize the library. Should be called from inside
            the application before any other call to the library. Should be called
            after the MPI is initialized.
  @param	filename        Name of the library configuration file.
  @param    module_name     Name of the module.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int farb_init(const char *filename, char *module_name)
{
  //  char* conf_filepath;
    int errno, mpi_initialized;
    char* s;
    int verbose;

    if(lib_initialized)
        return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        fprintf(stderr, "FARB Error: pFARB cannot be initialized before MPI is initialized. Aborting...\n");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &gl_my_rank);

    if(strlen(module_name)>MAX_COMP_NAME){
        fprintf(stderr, "FARB Error: module name %s too long\n", module_name);
        fflush(stderr);
        exit(1);
    }
    gl_conf.malloc_size = 0;

    gl_my_comp_name = (char*)farb_malloc(MAX_COMP_NAME);
    assert(gl_my_comp_name != NULL);
    strcpy(gl_my_comp_name, module_name);

    s = getenv("FARB_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_DBG_LEVEL;
    else
        gl_verbose = atoi(s);


    //during init only root will print out stuff
    if(gl_my_rank != 0){
        verbose = gl_verbose;
        gl_verbose = VERBOSE_ERROR_LEVEL;
    }
    /*Parse ini file and initialize components*/
    errno = load_config(filename, module_name);
    if(errno) goto panic_exit;

    /*Establish intercommunicators between components*/
    errno = init_comp_comm();
    if(errno) goto panic_exit;

    errno = init_data_distr();
    if(errno) goto panic_exit;

    errno = init_req_match_masters();
    if(errno) goto panic_exit;

    lib_initialized = 1;

    //enable print setting for other ranks again
    if(gl_my_rank != 0)
        gl_verbose = verbose;

    return 0;

panic_exit:

    farb_finalize();
    exit(1);
    return 1;
}

/**
  @brief	Function to finalize the library. Should be called from inside
            the application before the MPI is finalized.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int farb_finalize()
{
    int mpi_initialized;

    if(!lib_initialized) return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: pFARB cannot be finalized after MPI is finalized. Aborting...");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }

    finalize_comp_comm();

    clean_config();

    FARB_DBG(VERBOSE_DBG_LEVEL,"FARB: finalize");
    farb_free(gl_my_comp_name, MAX_COMP_NAME);
    if(gl_conf.masters != NULL)
        farb_free(gl_conf.masters, gl_conf.nmasters*sizeof(int));

    FARB_DBG(VERBOSE_DBG_LEVEL, "FARB memory leak size: %lu", gl_conf.malloc_size);
    lib_initialized = 0;
    fflush(stdout);
    fflush(stderr);
    return 0;


}


/*Interfaces to be used by a File I/O library*/
/**
  @brief	First checks if direct data transfer should be used for this file. If yes,
            writes a portion of data to corresponding memory buffer. If no, returns.
  @param	filename        file name for the memory buffer
  @param    offset          where in the file should this data be written
  @param    data_sz         size of the data to be written
  @param    data            pointer to the data to be written
  @return	number of bytes written

 */

_EXTERN_C_ void farb_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return;


    if(hdr_sz == 0){
        FARB_DBG(VERBOSE_DBG_LEVEL, "Header size for file %s is zero", filename);
        return;
    }
    write_hdr(fbuf, hdr_sz, header);
    return;
}

_EXTERN_C_ MPI_Offset farb_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;
    return read_hdr_chunk(fbuf, offset, chunk_sz, chunk);
}


_EXTERN_C_ void farb_create(const char *filename, int ncid)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_DBG_LEVEL, "File %s is not treated by FARB", filename);
        return;
    } else {
        FARB_DBG(VERBOSE_DBG_LEVEL, "Created file %s", filename);
    }
    fbuf->ncid = ncid;
}
/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void farb_open(const char *filename, MPI_Comm comm)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_DBG_LEVEL, "File %s is not treated by FARB", filename);
        return;
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "Opening file %s", filename);
    /*A hack for scale-letkf: if we are trying to open an alias,
    create an empty file, otherwise letkf will crash because
    it won't  find the file*/
    if( (fbuf->iomode == FARB_IO_MODE_MEMORY) && (strlen(fbuf->alias_name) > 0) && (strstr(filename, fbuf->alias_name) !=NULL)){
        MPI_File fh;
        int err;
        err = MPI_File_open(comm, (char*)filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
        CHECK_MPI(err);
        err = MPI_File_close(&fh);
        CHECK_MPI(err);
        FARB_DBG(VERBOSE_DBG_LEVEL, "Created a dummy alias file");

    }
    open_file(fbuf);
}

_EXTERN_C_ void farb_enddef(const char *filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return;

    if((fbuf->writer_id == gl_my_comp_id) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH))
        send_file_info(fbuf);
}
///**
//  @brief	Called before the corresponding file is opened. In case if components do normal File I/O
//            we need to synchronize the writer and reader(s) of the file so that reader(s) doesn't try to
//            open it before the writer finished writing.
//
//  @param	filename        file name for the memory buffer
//  @return	void
//
// */
//void farb_sync(char* filename)
//{
//
//}

/**
  @brief	Called when the corresponding file is closed. If the file is opened in write mode
            (output file) we will notify the reader that the file is ready to be transfered. If it's opened in read mode
            (input file), we will free the buffer.
            The mode is defined from the farb configuration file.
  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void farb_close(const char* filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->reader_id == gl_my_comp_id)
        MPI_Barrier(MPI_COMM_WORLD);
    close_file(fbuf);
}

/*called inside wait function in pnetcdf*/
_EXTERN_C_ int farb_match_ioreqs(const char* filename)
{
    int ret;
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;
    if(gl_conf.distr_mode != DISTR_MODE_REQ_MATCH) return 0;

    /*User will have to explicitly initiate matching*/
    if(fbuf->explicit_match) return 0;

    if(matching_flag)
        FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: farb_match_ioreqs is called for file %s, but a matching process has already started before.", fbuf->file_path);
    matching_flag = 1;
    ret = match_ioreqs(fbuf, 0);
    matching_flag = 0;

    return ret;
}

/*called by user to do explicit matching*/
/*
    User must specify either filename or ncid.
    intracomp_io_flag - if set to 1, matching of intracomponent io requests will be
    performed. This flag is intended for for situation when the writer component
    tries to read something from the file it is writing.
*/
_EXTERN_C_ int farb_match_io(const char *filename, int ncid, int intracomp_io_flag )//, int match_all)
{
    if(!lib_initialized) return 0;
    if(gl_conf.distr_mode != DISTR_MODE_REQ_MATCH) return 0;
//    if(match_all){
//        file_buffer_t *fbuf = gl_filebuf_list;
//        while(fbuf != NULL){
//            if(fbuf->iomode == FARB_IO_MODE_MEMORY)
//                match_ioreqs(fbuf);
//            fbuf = fbuf->next;
//        }
//
//    } else{
        file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
        if(fbuf == NULL){

            if( (filename != NULL) && (strlen(filename) == 0) )
                FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: file with ncid %d is not treated by FARB (not in configuration file). Explicit matching ignored.", ncid);
            else
                FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: file %s (ncid %d) is not treated by FARB (not in configuration file). Explicit matching ignored.", filename, ncid);
            return 0;
        } else {
            FARB_DBG(VERBOSE_DBG_LEVEL, "file %s, ncid %d, io flag %d", fbuf->file_path, ncid, intracomp_io_flag);
        }
        if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;
        if(!fbuf->explicit_match){
            FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: calling farb_match_io but explicit match for file %s not enabled. Ignored.", filename);
            return 0;
        }

        if( intracomp_io_flag && (gl_my_comp_id != fbuf->writer_id) ){
            FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: farb_match_io: intracomp_io_flag(%d) can only be set for the writer component", intracomp_io_flag);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        if(matching_flag){
            FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: farb_match_io is called for file %s, but a matching process has already started before.", fbuf->file_path);
            if(intracomp_io_flag){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: farb_match_io: intracomp_io_flag is set but the process is already matching io. This is not allowed.");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
        }
        matching_flag = 1;
        match_ioreqs(fbuf, intracomp_io_flag);
        matching_flag = 0;
  // }
    return 0;
}
/* The user has to state whether the process needs to match all read or write requests.
   Because the process of matching for a reader and writer is not the same. */
_EXTERN_C_ void farb_match_io_all(int rw_flag)
{
    if(!lib_initialized) return;
    if(gl_conf.distr_mode != DISTR_MODE_REQ_MATCH) return;

    if(rw_flag == FARB_READ){
        FARB_DBG(VERBOSE_WARNING_LEVEL, "farb_match_io_all() cannot be used in processes that read files. Ignoring.");
        return;
    }
    if(matching_flag)
        FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: farb_match_io_all is called, but a matching process has already started before.");
    matching_flag = 1;
    match_ioreqs_all(rw_flag);
    matching_flag = 0;
    return;
}

/**
    @brief  Check if the file is intended to be written by this component
    @param  filename    name of the file
    @return 1 - yes, 0 - no
*/
_EXTERN_C_ int farb_write_flag(const char* filename)
{
    if(!lib_initialized)return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    return (fbuf->writer_id == gl_my_comp_id) ? 1 : 0;
}

/**
    @brief  Check if the file intended to be read by this component
    @param  filename    name of the file
    @return 1 - yes, 0 - no
*/
_EXTERN_C_ int farb_read_flag(const char* filename)
{
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    return (fbuf->reader_id == gl_my_comp_id) ? 1 : 0;
}

_EXTERN_C_ MPI_Offset farb_read_write_var(const char *filename,
                                          int varid,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          const MPI_Offset *stride,
                                          const MPI_Offset *imap,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag,
                                          int *request)
{
    MPI_Offset ret;

    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;

    if(boundary_check(fbuf, varid, start, count ))
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);

    if( rw_flag != FARB_READ && rw_flag != FARB_WRITE){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "rw_flag value incorrect (%d)", rw_flag);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(rw_flag==FARB_WRITE && fbuf->reader_id == gl_my_comp_id){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: reader component cannot write to the file %s", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    switch(gl_conf.distr_mode){
            case DISTR_MODE_STATIC:
                ret = buf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag);
                break;
            case DISTR_MODE_REQ_MATCH:
                if(request == NULL)
                    ret = nbuf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag, NULL);
                else
                    ret = nbuf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag, request);
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: 2 unknown data distribution mode");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    return ret;
}


/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/
//TODO do we really need it?
_EXTERN_C_ int farb_io_mode(const char* filename)
{
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL)
        return FARB_IO_MODE_UNDEFINED;
    return fbuf->iomode;
}

_EXTERN_C_ int farb_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i, ret;
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;
//TODO make this work for fixed distr matching as well
    /*Check if there is an UNLIMITED dimension: currently it's not supported*/
    for(i = 0; i < ndims; i++)
        FARB_DBG(VERBOSE_DBG_LEVEL, "varid %d, dim %d size %llu", varid, i, shape[i]);
    /*For now, can only support unlimited dimension if it's the first dimension array*/
    if( (ndims > 0) && (shape[0] == FARB_UNLIMITED))
        FARB_DBG(VERBOSE_DBG_LEVEL, "var has unlimited dimension");
    for(i = 1; i < ndims; i++){
        //we can support unlimited dimension if it's the first dimension
        //else abort
        if(shape[i] == FARB_UNLIMITED){
            FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: currently cannot support when unlimited dimension is not in the slowest changing dimension (dim 0). Aborting.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(dtype == MPI_DATATYPE_NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: datatype for var %d (file %s) is null. Aborting.", varid, fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(frt_indexing && ndims > 0){
        MPI_Offset *cshape = (MPI_Offset*)farb_malloc(sizeof(MPI_Offset)*ndims);
        assert(cshape != NULL);

        for(i = 0; i < ndims; i++)
            if(shape[i] != FARB_UNLIMITED)
                cshape[i] = shape[i] + 1;
            else
                cshape[i] = 0;
        ret = def_var(fbuf, varid, ndims, dtype, cshape);
        farb_free(cshape, sizeof(MPI_Offset)*ndims);
    } else
      ret = def_var(fbuf, varid, ndims, dtype, shape);

    return ret;
}

_EXTERN_C_ int farb_set_distr_count(const char* filename, int varid, int count[])
{
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != FARB_IO_MODE_MEMORY) return 0;
    return set_distr_count(fbuf, varid, count);
}

/*  Fortran Interfaces  */

void farb_init_(const char *filename, char *module_name, int* ierr)
{
    frt_indexing = 1;
    FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: Dealing with fortran. Increment each dimension by 1");
    *ierr = farb_init(filename, module_name);
}

void farb_finalize_(int* ierr)
{
    *ierr = farb_finalize();
}

void farb_match_io_(const char *filename, int *ncid, int *intracomp_io_flag, int *ierr)
{
    *ierr = farb_match_io(filename, *ncid, *intracomp_io_flag);
}

void farb_match_io_all_(int *rw_flag)
{
    farb_match_io_all(*rw_flag);
}
