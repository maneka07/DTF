#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>


//TODO figure out with hierarchical way of including headers
#include "pfarb_common.h"
#include "pfarb.h"
#include "pfarb_init_finalize.h"
#include "pfarb_util.h"
#include "pfarb_buf_io.h"


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
//_EXTERN_C_ size_t farb_write(const char* filename, off_t const offset, const size_t data_sz, void *const data)
//{
////    if(!lib_initialized) return 0;
////    if(farb_io_mode(filename) != FARB_IO_MODE_MEMORY) return 0;
////    return mem_write(filename, offset, data_sz, data);

    //return 0;
//}

///**
  //@brief	First checks if direct data transfer should be used for this file. If yes,
            //reads a portion of data to corresponding memory buffer. If no, returns.
  //@param	filename        file name for the memory buffer
  //@param    offset          offset at which to read the data
  //@param    data_sz         size of the data to be read
  //@param    data            pointer to the buffer where to read the data to
  //@return	number of bytes read

 //*/
//_EXTERN_C_ size_t farb_read(const char *filename, off_t const offset, const size_t data_sz, void *const data)
//{

////    if(!lib_initialized) return 0;
////    if(farb_io_mode(filename) != FARB_IO_MODE_MEMORY) return 0;
////    FARB_DBG(VERBOSE_DBG_LEVEL,   "read %s", filename);
////    return mem_read(filename, offset, data_sz, data);

    //return 0;
//}


_EXTERN_C_ void farb_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    if(!lib_initialized) return;
    //if(farb_io_mode(filename) != FARB_IO_MODE_MEMORY) return;
    if(hdr_sz == 0){
        FARB_DBG(VERBOSE_DBG_LEVEL, "Header size for file %s is zero", filename);
        return;
    }
    write_hdr(filename, hdr_sz, header);
    return;
}

_EXTERN_C_ MPI_Offset farb_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
     if(!lib_initialized) return 0;

     return read_hdr_chunk(filename, offset, chunk_sz, chunk);
}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void farb_open(const char *filename)
{
    if(!lib_initialized) return;
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Enter farb_open %s", filename);
    if(get_read_flag(filename)){
        while(!file_buffer_ready(filename))
            //force progress until the writer finishes with the file.
            progress_io();
//        switch(gl_conf.distr_mode){
//            case DISTR_MODE_STATIC:
//                while(!file_buffer_ready(filename))
//                    //force progress until the writer finishes with the file.
//                    progress_io();
//                break;
//            case DISTR_MODE_NONBUFFERED_REQ_MATCH:
//
//                break;
//            default:
//                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: unknown data distribution mode");
//                MPI_Abort(0);
//        }
    }

    else
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Exit farb_open %s", filename);
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
    FARB_DBG(VERBOSE_DBG_LEVEL,   "close %s", filename);
    close_file(filename);
}

_EXTERN_C_ int farb_match_ioreqs(const char* filename)
{
    if(!lib_initialized) return 0;
    return match_ioreqs(filename);
}

/**
    @brief  Check if the file is intended to be written by this component
    @param  filename    name of the file
    @return 1 - yes, 0 - no
*/
_EXTERN_C_ int farb_write_flag(const char* filename)
{
    int flag;
    if(!lib_initialized)return 0;

    flag = get_write_flag(filename);
    return flag;
}

/**
    @brief  Check if the file intended to be read by this component
    @param  filename    name of the file
    @return 1 - yes, 0 - no
*/
_EXTERN_C_ int farb_read_flag(const char* filename)
{
    if(!lib_initialized) return 0;
    return get_read_flag(filename);
}

_EXTERN_C_ MPI_Offset farb_read_write_var(const char *filename, int varid, const MPI_Offset *start, const MPI_Offset *count,
                                          const MPI_Offset *stride, const MPI_Offset *imap, MPI_Datatype dtype, void *buf, int rw_flag)
{
    if(!lib_initialized) return 0;

    return read_write_var(filename, varid, start, count, stride, imap, dtype, buf, rw_flag);
}

/**
    @brief  Progress with sending/receiving of data
    @return void
*/
_EXTERN_C_ void farb_progress_io()
{

    if(!lib_initialized) return;
    progress_io();
}

/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/
_EXTERN_C_ int farb_io_mode(const char* filename)
{
    if(!lib_initialized) return 0;
    return get_io_mode(filename);
}

_EXTERN_C_ int farb_def_var(const char* filename, int varid, int ndims, MPI_Offset *shape)
{
    if(!lib_initialized) return 0;

   // if(farb_io_mode(filename) != FARB_IO_MODE_MEMORY) return 0;

    return def_var(filename, varid, ndims, shape);
}

_EXTERN_C_ int farb_set_distr_count(const char* filename, int varid, int count[])
{
    if(!lib_initialized) return 0;

    return set_distr_count(filename, varid, count);
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
