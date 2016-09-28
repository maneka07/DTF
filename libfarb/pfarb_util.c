/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */


#include <assert.h>
#include <string.h>
#include <mpi.h>
#include "pfarb_util.h"
#include "pfarb_mem.h"
#include "pfarb.h"
#include "pfarb_buf_io.h"

int get_write_flag(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->write_flag;
}

int get_read_flag(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }

    return fbuf->read_flag;
}

int file_buffer_ready(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->is_ready;
}


void progress_io()
{
    MPI_Status status;
    int i, flag, src, errno;
    char filename[MAX_FILE_NAME];
    file_buffer_t *fbuf;

    for(i = 0; i < gl_ncomp; i++){
        if( i == gl_my_comp_id || gl_comps[i].intercomm == MPI_COMM_NULL)
            continue;

        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[i].intercomm, &flag, &status);
        if(!flag)
            continue;

        switch(status.MPI_TAG){
            case FILE_READY_TAG:
                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[i].intercomm, &status);
                CHECK_MPI(errno);
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);

                fbuf = find_file_buffer(gl_filebuf_list, filename);
                assert(fbuf != NULL);

                if(fbuf->mode == FARB_IO_MODE_MEMORY){
                    /*Notify I am ready to receive*/
                    errno = MPI_Send(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm);
                    CHECK_MPI(errno);
                    receive_data(fbuf, src, gl_comps[i].intercomm);
                    fbuf->distr_ndone++;
                } else
                    fbuf->distr_ndone++;
                if(fbuf->distr_ndone == fbuf->distr_nranks)
                    fbuf->is_ready = 1;
                break;
            case RECV_READY_TAG:

                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm, &status);

                CHECK_MPI(errno);
                fbuf = find_file_buffer(gl_filebuf_list, filename);
                assert(fbuf != NULL);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Receive RECV_READY_TAG notif for %s", filename);

                send_data(fbuf, src, gl_comps[i].intercomm);

                fbuf->distr_ndone++;
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown tag %d", status.MPI_TAG);
                assert(0);
        }
    }
}

MPI_Offset to_1d_index(int ndims, const MPI_Offset *shape, MPI_Offset *coord)
{
      int i, j;
      MPI_Offset idx = 0, mem=0;

      if(ndims == 0) //scalar
        return 0;
      else if(ndims == 1) //1d array
        return *coord;

      for(i = 0; i < ndims; i++){
        mem = coord[i];
        for(j = i+1; j < ndims; j++)
          mem *= shape[j];
        idx += mem;
      }
    return idx;
}

void notify_file_ready(const char* filename)
{
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);
    int comp_id, i, dest, errno;

    comp_id = fbuf->reader_id;
    for(i = 0; i < fbuf->distr_nranks; i++){
        dest = fbuf->distr_ranks[i];
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader rank %d that file %s is ready", dest, fbuf->file_path);
        errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, dest, FILE_READY_TAG, gl_comps[comp_id].intercomm);
        CHECK_MPI(errno);

        if(fbuf->mode == FARB_IO_MODE_FILE)
            fbuf->distr_ndone++;
    }

}

void close_file(const char *filename)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->write_flag){
        notify_file_ready(filename);
        while(fbuf->distr_ndone != fbuf->distr_nranks)
            progress_io();
    }

    delete_file_buffer(&gl_filebuf_list, fbuf);
}

int get_io_mode(const char* filename)
{

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL)
        return FARB_IO_MODE_UNDEFINED;
    return fbuf->mode;
}

int def_var(const char* filename, int varid, int ndims, MPI_Offset *shape)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    farb_var_t *var = new_var(varid, ndims, shape);

    add_var(&buf->vars, var);
    buf->var_cnt++;
    return 0;
}

/*Write pnetcdf header*/
void write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    buf->hdr_sz = hdr_sz;
    buf->header = malloc(hdr_sz);
    assert(buf->header != NULL);
    memcpy(buf->header, header, (size_t)hdr_sz);
    return;
}

/*Read pnetcdf header*/
MPI_Offset read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    if(offset+chunk_sz > buf->hdr_sz){
        FARB_DBG(VERBOSE_DBG_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, buf->hdr_sz);
        chunk_sz = buf->hdr_sz - offset;
    }

    memcpy(chunk, buf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}

/*Set how many elements in each dimension to distribute to ranks.
  File's distrib_pattern must be set to scatter */
int set_distr_count(const char* filename, int varid, int count[])
{
    int i;
    MPI_Offset *cnt;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->mode != FARB_IO_MODE_MEMORY){
        //nothing to do
        return 0;
    }

    if(fbuf->distr_pattern != DISTR_PATTERN_SCATTER){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: cannot set distribute count. File's distribute pattern must be <scatter>.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: No variable with such id (%d)", varid);
        assert(0);
    }

    cnt = (MPI_Offset*)realloc(var->distr_count, var->ndims*sizeof(MPI_Offset));
    assert(cnt != NULL);
    var->distr_count = cnt;
    for(i = 0; i < var->ndims; i++)
        var->distr_count[i] = (MPI_Offset)count[i];

    return 0;
}
