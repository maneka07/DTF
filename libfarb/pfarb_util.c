/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */


#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "pfarb_util.h"
#include "pfarb_mem.h"
#include "pfarb.h"
#include "pfarb_buf_io.h"
#include "pfarb_file_buffer.h"
#include "pfarb_common.h"
#include "pfarb_req_match.h"

//int get_write_flag(const char* filename)
//{
//    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//    if(fbuf == NULL){
//        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
//        assert(0);
//    }
//    if(fbuf->writer_id == gl_my_comp_id)
//        return 1;
//    else
//        return 0;
//}
//

//int get_read_flag(const char* filename)
//{
//    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//    if(fbuf == NULL){
//        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
//        assert(0);
//    }
//
//    if(fbuf->reader_id == gl_my_comp_id)
//        return 1;
//    else
//        return 0;
//}

int file_buffer_ready(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->is_ready;
}


/*Used in static data distribution*/
void progress_io()
{
    MPI_Status status;
    int i, flag, src, errno;
    char filename[MAX_FILE_NAME];
    file_buffer_t *fbuf;
    int bufsz;

    for(i = 0; i < gl_ncomp; i++){
        if( gl_comps[i].comm == MPI_COMM_NULL)
            continue;

        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[i].comm, &flag, &status);
        if(!flag)
            continue;

        switch(status.MPI_TAG){
            case FILE_READY_TAG:
                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[i].comm, &status);
                CHECK_MPI(errno);
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);

                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                assert(fbuf != NULL);

                if(fbuf->iomode == FARB_IO_MODE_MEMORY){
                    /*Notify I am ready to receive*/
                    errno = MPI_Send(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].comm);
                    CHECK_MPI(errno);

                    /*Receive the header*/
                    MPI_Probe(src, HEADER_TAG, gl_comps[i].comm, &status);
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    fbuf->hdr_sz = (MPI_Offset)bufsz;
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Hdr size to receive %d", bufsz);
                    fbuf->header = malloc((size_t)bufsz);
                    assert(fbuf->header != NULL);
                    errno = MPI_Recv(fbuf->header, bufsz, MPI_BYTE, src, HEADER_TAG, gl_comps[i].comm, &status);
                    CHECK_MPI(errno);

                    receive_data(fbuf, src, gl_comps[i].comm);
                    fbuf->distr_ndone++;
                    if(fbuf->distr_ndone == fbuf->distr_nranks)
                        fbuf->is_ready = 1;
                } else
                    fbuf->is_ready = 1;

                break;
            case RECV_READY_TAG:

                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].comm, &status);

                CHECK_MPI(errno);
                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                assert(fbuf != NULL);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Receive RECV_READY_TAG notif for %s", filename);
                /*Send the hearder*/
                errno = MPI_Send(fbuf->header, (int)fbuf->hdr_sz, MPI_BYTE, src, HEADER_TAG, gl_comps[i].comm);
                CHECK_MPI(errno);
                send_data(fbuf, src, gl_comps[i].comm);

                fbuf->distr_ndone++;
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown tag %d", status.MPI_TAG);
                assert(0);
        }
    }
}

MPI_Offset last_1d_index(int ndims, const MPI_Offset *shape)
{
    MPI_Offset ret = 0;

    if(ndims > 0){
        int i;
        MPI_Offset *tmp = (MPI_Offset*)malloc(ndims*sizeof(MPI_Offset));
        assert(tmp != NULL);
        for(i = 0; i < ndims; i++)
            tmp[i] = shape[i] - 1;
        ret = to_1d_index(ndims, shape, tmp);
        free(tmp);
    }

    return ret;
}

MPI_Offset to_1d_index(int ndims, const MPI_Offset *shape, const MPI_Offset *coord)
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

void notify_file_ready(file_buffer_t *fbuf)
{
    int comp_id, i, dest, errno, nranks;

    comp_id = fbuf->reader_id;
    if(fbuf->iomode == FARB_IO_MODE_FILE){
            MPI_Comm_remote_size(gl_comps[comp_id].comm, &nranks);
            if(gl_my_rank < nranks){
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader rank %d that file %s is ready", gl_my_rank, fbuf->file_path);
                errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, gl_my_rank, FILE_READY_TAG, gl_comps[comp_id].comm);
                CHECK_MPI(errno);
            }
    } else {
        switch(gl_conf.distr_mode){
            case DISTR_MODE_STATIC:
                for(i = 0; i < fbuf->distr_nranks; i++){
                    dest = fbuf->distr_ranks[i];
                    FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader rank %d that file %s is ready", dest, fbuf->file_path);
                    errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, dest, FILE_READY_TAG, gl_comps[comp_id].comm);
                    CHECK_MPI(errno);

                    if(fbuf->iomode == FARB_IO_MODE_FILE)
                        fbuf->distr_ndone++;
                }
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: 4 unknown data distribution mode");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }


}

void close_file(file_buffer_t *fbuf)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", fbuf->file_path);
//TODO for every else for writer_id write explicit compare of reader id (in case have 3+ components)
    if(fbuf->writer_id == gl_my_comp_id){

        if((fbuf->iomode == FARB_IO_MODE_FILE) && (gl_conf.distr_mode == DISTR_MODE_STATIC)){
            notify_file_ready(fbuf);
            while(fbuf->distr_ndone != fbuf->distr_nranks)
                progress_io();
        } else if(fbuf->iomode == FARB_IO_MODE_FILE) {
            notify_file_ready(fbuf);
        } else if((fbuf->iomode == FARB_IO_MODE_MEMORY) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH)){
            while( !fbuf->fclosed_flag)
                progress_io_matching();
            assert(fbuf->ioreq_cnt == 0);
        }
    } else if (fbuf->reader_id == gl_my_comp_id){
        if((fbuf->iomode == FARB_IO_MODE_MEMORY) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH)){
            /*Reader ranks synch to ensure that all read requests have been completed*/
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Synchronize before closing the file %s", fbuf->file_path);
//            MPI_Barrier(gl_comps[gl_my_comp_id].comm);
            assert(fbuf->ioreq_cnt == 0);
            /*Notify writers they can delete their write requests*/
            errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, gl_conf.my_master, IO_CLOSE_FILE_TAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(errno);
        }
    }

    delete_file_buffer(&gl_filebuf_list, fbuf);
}

void open_file(file_buffer_t *fbuf)
{
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Enter farb_open %s", fbuf->file_path);
    if(fbuf->reader_id == gl_my_comp_id){
        while(!fbuf->is_ready)
            if(fbuf->iomode == FARB_IO_MODE_FILE)
                progress_io();
            else if(fbuf->iomode == FARB_IO_MODE_MEMORY){
                //force progress until the writer finishes with the file.
                switch(gl_conf.distr_mode){
                    case DISTR_MODE_STATIC:
                        progress_io();
                        break;
                    case DISTR_MODE_REQ_MATCH:
                        progress_io_matching();
                        break;
                    default:
                        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: 1 unknown data distribution mode");
                        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
                }
            }
    }
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Exit farb_open %s", fbuf->file_path);
}

int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Offset el_sz, MPI_Offset *shape)
{
    farb_var_t *var = new_var(varid, ndims, el_sz, shape);
    add_var(&fbuf->vars, var);
    fbuf->var_cnt++;
    return 0;
}

/*Write pnetcdf header*/
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Writing header (sz %d)", (int)hdr_sz);
    fbuf->hdr_sz = hdr_sz;
    fbuf->header = malloc(hdr_sz);
    assert(fbuf->header != NULL);
    memcpy(fbuf->header, header, (size_t)hdr_sz);
    return;
}

/*Read pnetcdf header*/
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(offset+chunk_sz > fbuf->hdr_sz){
        FARB_DBG(VERBOSE_ALL_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, fbuf->hdr_sz);
        chunk_sz = fbuf->hdr_sz - offset;
    }

    memcpy(chunk, (unsigned char*)fbuf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}
