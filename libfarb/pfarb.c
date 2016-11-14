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
//TODO change when request is sent to a master
//TODO unlimited dimension!!
//1) how do we know current size?

int lib_initialized=0;
int gl_verbose;
int gl_my_rank;
struct farb_config gl_conf;
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

    gl_my_comp_name = (char*)malloc(MAX_COMP_NAME);
    assert(gl_my_comp_name != NULL);
    strcpy(gl_my_comp_name, module_name);

    s = getenv("FARB_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_ERROR_LEVEL;
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
    free(gl_my_comp_name);
    if(gl_conf.masters != NULL)
        free(gl_conf.masters);

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
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return;


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
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;
    return read_hdr_chunk(fbuf, offset, chunk_sz, chunk);
}


_EXTERN_C_ void farb_create(const char *filename, int ncid)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ALL_LEVEL, "File %s is not treated by FARB", filename);
        return;
    }
    fbuf->ncid = ncid;
}
/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void farb_open(const char *filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ALL_LEVEL, "File %s is not treated by FARB", filename);
        return;
    }
    open_file(fbuf);
}

_EXTERN_C_ void farb_enddef(const char *filename)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return;

    if((fbuf->writer_id == gl_my_comp_id) && (gl_conf.distr_mode == DISTR_MODE_NONBUFFERED_REQ_MATCH))
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
    close_file(fbuf);
}

/*called inside wait function in pnetcdf*/
_EXTERN_C_ int farb_match_ioreqs(const char* filename)
{
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;
    if(gl_conf.distr_mode != DISTR_MODE_NONBUFFERED_REQ_MATCH) return 0;

    /*User will have to explicitly initiate matching*/
    if(gl_conf.explicit_match) return 0;

    return match_ioreqs(fbuf);
}

/*called by user to do explicit matching*/
_EXTERN_C_ int farb_match_io(const char *filename, int match_all)
{
    if(!lib_initialized) return 0;
    if(gl_conf.distr_mode != DISTR_MODE_NONBUFFERED_REQ_MATCH) return 0;
    if(!gl_conf.explicit_match) return 0;
    if(match_all){
        file_buffer_t *fbuf = gl_filebuf_list;
        while(fbuf != NULL){
            if(fbuf->mode == FARB_IO_MODE_MEMORY)
                match_ioreqs(fbuf);
            fbuf = fbuf->next;
        }

    } else{
        file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
        if(fbuf == NULL) return 0;
        if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;
        match_ioreqs(fbuf);
    }
    return 0;
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
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;

    switch(gl_conf.distr_mode){
            case DISTR_MODE_STATIC:
                ret = buf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag);
                break;
            case DISTR_MODE_NONBUFFERED_REQ_MATCH:
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
    return fbuf->mode;
}

_EXTERN_C_ int farb_def_var(const char* filename, int varid, int ndims, MPI_Offset el_sz, MPI_Offset *shape)
{
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;
    FARB_DBG(VERBOSE_ALL_LEVEL, "el_sz %d", (int)el_sz);
    return def_var(fbuf, varid, ndims, el_sz, shape);
}

_EXTERN_C_ int farb_set_distr_count(const char* filename, int varid, int count[])
{
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->mode != FARB_IO_MODE_MEMORY) return 0;
    return set_distr_count(fbuf, varid, count);
}

/*  Fortran Interfaces  */

void farb_init_(const char *filename, char *module_name, int* ierr)
{
    *ierr = farb_init(filename, module_name);
}

void farb_finalize_(int* ierr)
{
    *ierr = farb_finalize();
}

void farb_match_io_(const char *filename, int match_all, int *ierr)
{
    *ierr = farb_match_io(filename, match_all);
}
