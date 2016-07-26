#include <mpi.h>
#include <stdio.h>

#include "farb.h"
#include "farb_util_new.h"
#include <assert.h>

//TODO change everithing file size related to long long and check casts

int lib_initialized=0;
int gl_verbose;
int gl_my_rank;
farb_settings_t gl_sett;
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
    MPI_Datatype dt[2] = {MPI_CHAR, MPI_UNSIGNED};
    int blocklen[2] = {MAX_FILE_NAME, 1};
    MPI_Aint displ[2];

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

    strcpy(gl_my_comp_name, module_name);

    s = getenv("FARB_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_ERROR_LEVEL;
    else
        gl_verbose = atoi(s);
    //for now only root will print out stuff
    //if(gl_my_rank != 0)
      //  gl_verbose = VERBOSE_NONE;

    s = getenv("FARB_NODE_SZ");
    if(s == NULL){
        gl_sett.node_sz = DEFAULT_BUFFER_NODE_SIZE;
    } else
        gl_sett.node_sz = atoi(s)*1024;

    FARB_DBG(VERBOSE_DBG_LEVEL,   "Buffer node size set to %ld bytes", gl_sett.node_sz);

    s = getenv("FARB_MESSAGE_SZ");
    if(s == NULL){
        gl_sett.msg_sz = (int)gl_sett.node_sz;
    } else
        gl_sett.node_sz = atoi(s)*1024;

    if(gl_sett.msg_sz > (int) gl_sett.node_sz){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: message size (%d) for direct data transfer cannot exceed the memory node size. Setting to same as the node size (%d)", gl_sett.msg_sz, (int)gl_sett.node_sz);
        gl_sett.msg_sz = (int)gl_sett.node_sz;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL,   "MPI message size set to %d bytes", gl_sett.msg_sz);

    /*Parse ini file and initialize components*/
    errno = load_config(filename, module_name);
    if(errno) goto panic_exit;

    /*Establish intercommunicators between components*/
    errno = init_comp_comm();
    if(errno) goto panic_exit;


    /*Create MPI_datatype for file-related notification messages*/
    displ[0] = offsetof(msg_ready_notif_t, filename);
    displ[1] = offsetof(msg_ready_notif_t, file_sz);

    errno = MPI_Type_create_struct(2, blocklen, displ, dt, &msg_ready_datatype);
    assert(errno == MPI_SUCCESS);
    errno = MPI_Type_commit(&msg_ready_datatype);
    assert(errno == MPI_SUCCESS);

    lib_initialized = 1;
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

    MPI_Type_free(&msg_ready_datatype);
    finalize_comp_comm();

    clean_config();

    FARB_DBG(VERBOSE_DBG_LEVEL,   "Farb: finalize\n");
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
_EXTERN_C_ size_t farb_write(const char* filename, off_t const offset, const size_t data_sz, void *const data)
{
    if(!lib_initialized) return 0;
    return mem_write(filename, offset, data_sz, data);
}

/**
  @brief	First checks if direct data transfer should be used for this file. If yes,
            reads a portion of data to corresponding memory buffer. If no, returns.
  @param	filename        file name for the memory buffer
  @param    offset          offset at which to read the data
  @param    data_sz         size of the data to be read
  @param    data            pointer to the buffer where to read the data to
  @return	number of bytes read

 */
_EXTERN_C_ size_t farb_read(const char* filename, off_t const offset, const size_t data_sz, void *const data)
{

    if(!lib_initialized) return 0;
    FARB_DBG(VERBOSE_DBG_LEVEL,   "read %s", filename);
    return mem_read(filename, offset, data_sz, data);

}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void farb_open(const char* filename)
{
    if(!lib_initialized) return;
    FARB_DBG(VERBOSE_DBG_LEVEL,   "open %s", filename);

    if(get_read_flag(filename))
        while(!file_buffer_ready(filename)){
            //force progress until the writer finishes with the file.
            progress_io();
        }
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
    FARB_DBG(VERBOSE_DBG_LEVEL,   "write flag %s is %d", filename, flag);
    return flag;
}

/**
    @brief  Check if the file intended to be read by this component
    @param  filename    name of the file
    @return 1 - yes, 0 - no
*/
_EXTERN_C_ int farb_read_flag(const char* filename)
{
    int flag;

    if(!lib_initialized) return 0;

    flag =  get_read_flag(filename);
    FARB_DBG(VERBOSE_DBG_LEVEL,   "read flag %s is %d", filename, flag);
    return flag;
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
    return get_io_mode(filename);
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
