/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */


#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "dtf_util.h"
#include "dtf_mem.h"
#include "dtf.h"
#include "dtf_file_buffer.h"
#include "dtf_common.h"
#include "dtf_req_match.h"

extern file_info_req_q_t *gl_finfo_req_q;

//int get_write_flag(const char* filename)
//{
//    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//    if(fbuf == NULL){
//        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: file %s not in the configuration file", filename);
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
//        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: file %s not in the configuration file", filename);
//        assert(0);
//    }
//
//    if(fbuf->reader_id == gl_my_comp_id)
//        return 1;
//    else
//        return 0;
//}



/*global rank 0 periodically checks if it can
  process a file info request from a reader*/
void process_file_info_req_queue()
{
    if(gl_finfo_req_q == NULL)
        return;

    file_buffer_t *fbuf;
    file_info_req_q_t *req = gl_finfo_req_q, *tmp;

    while(req != NULL){
        fbuf = find_file_buffer(gl_filebuf_list, req->filename, -1);
        assert(fbuf != NULL);
        if(fbuf->root_writer == -1){
            req = req->next;
            continue;
        }

        if(gl_my_rank == fbuf->root_writer){
            DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info req for file %s", req->filename);
            memcpy(&fbuf->root_reader, (unsigned char*)(req->buf)+MAX_FILE_NAME, sizeof(int));

            if(fbuf->iomode == DTF_IO_MODE_MEMORY && gl_conf.distr_mode == DISTR_MODE_REQ_MATCH){
                memcpy(&(fbuf->mst_info->nrranks_opened), (unsigned char*)(req->buf)+MAX_FILE_NAME+sizeof(int), sizeof(int));
                send_file_info(fbuf, fbuf->root_reader);
            } else if(fbuf->iomode == DTF_IO_MODE_FILE){
                DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");
                errno = MPI_Send(&fbuf->ncid,1, MPI_INT, fbuf->root_reader, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm);
                CHECK_MPI(errno);
            }
        } else {
            int errno;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Forward the request to root writer %d", fbuf->root_writer);
            /*Forward this request to the root master*/
            errno = MPI_Send(req->buf, (int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm);
            CHECK_MPI(errno);
        }

        tmp = req->next;
        //dequeue
        if(req == gl_finfo_req_q){
            gl_finfo_req_q = tmp;
        } else{
            req->prev->next = req->next;
            if(req->next != NULL)
                req->next->prev = req->prev;
        }
        dtf_free(req->buf, MAX_FILE_NAME+sizeof(int)*2);
        dtf_free(req, sizeof(file_info_req_q_t));
        req = tmp;
        DTF_DBG(VERBOSE_ALL_LEVEL, "Freed request");
    }
}

int file_buffer_ready(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->is_ready;
}


//Used in static data distribution
//void progress_io()
//{
//    MPI_Status status;
//    int i, flag, src, errno;
//    char filename[MAX_FILE_NAME];
//    file_buffer_t *fbuf;
//    int bufsz;
//    void *rbuf;
//
//    /*first check if there are any requests for file info
//      that we can process. */
//    process_file_info_req_queue();
//
//    for(i = 0; i < gl_ncomp; i++){
//        if( gl_comps[i].comm == MPI_COMM_NULL)
//            continue;
//
//        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[i].comm, &flag, &status);
//        if(!flag)
//            continue;
//
//        src = status.MPI_SOURCE;
//
//        switch(status.MPI_TAG){
//            case FILE_INFO_REQ_TAG:
//                DTF_DBG(VERBOSE_DBG_LEVEL, "Recv root req from rank %d (comp %s)", src, gl_comps[i].name);
//                rbuf = dtf_malloc(MAX_FILE_NAME+2*sizeof(int));
//                assert(rbuf != NULL);
//                errno = MPI_Recv(rbuf, (int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, src, FILE_INFO_REQ_TAG, gl_comps[i].comm, &status);
//                CHECK_MPI(errno);
//                memcpy(filename, rbuf, MAX_FILE_NAME);
//                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//                assert(fbuf != NULL);
//
//                /* If the global master zero (global rank 0) couldn't process the
//                   request immediately because it doesn't know yet who is the root
//                   writer for the file, the req has to be enqueue and, later, the master
//                   will periodically check the queue to see if it can process the
//                   request.
//                 */
//                if(fbuf->root_writer == -1){
//                    DTF_DBG(VERBOSE_DBG_LEVEL, "Don't know who is the root for %s now. Will queue the req.", fbuf->file_path);
//                    file_info_req_q_t *req = dtf_malloc(sizeof(file_info_req_q_t));
//                    assert(req != NULL);
//                    memcpy(req->filename, filename, MAX_FILE_NAME);
//                    req->buf = rbuf;
//                    req->next = NULL;
//                    req->prev = NULL;
//
//                    enqueue
//                    if(gl_finfo_req_q == NULL)
//                        gl_finfo_req_q = req;
//                    else {
//                        req->next = gl_finfo_req_q;
//                        gl_finfo_req_q->prev = req;
//                        gl_finfo_req_q = req;
//                    }
//                    break;
//                }
//                if(gl_my_rank == fbuf->root_writer){
//                    DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");
//                    memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
//
//                    errno = MPI_Send(&fbuf->ncid,1, MPI_INT, fbuf->root_reader, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm);
//                    CHECK_MPI(errno);
//
//                } else {
//                    DTF_DBG(VERBOSE_DBG_LEVEL, "Forward the request to rank %d", fbuf->root_writer);
//                    Forward this request to the root master
//                    errno = MPI_Send(rbuf,(int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm);
//                    CHECK_MPI(errno);
//                }
//                dtf_free(rbuf, MAX_FILE_NAME+2*sizeof(int));
//                break;
//            case FILE_READY_TAG:
//                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[i].comm, &status);
//                CHECK_MPI(errno);
//                DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);
//
//                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//                assert(fbuf != NULL);
//
//                if(fbuf->iomode == DTF_IO_MODE_MEMORY){
//                    Notify I am ready to receive
//                    errno = MPI_Send(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].comm);
//                    CHECK_MPI(errno);
//
//                    Receive the header
//                    MPI_Probe(src, HEADER_TAG, gl_comps[i].comm, &status);
//                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
//                    fbuf->hdr_sz = (MPI_Offset)bufsz;
//                    DTF_DBG(VERBOSE_DBG_LEVEL, "Hdr size to receive %d", bufsz);
//                    fbuf->header = dtf_malloc((size_t)bufsz);
//                    assert(fbuf->header != NULL);
//                    errno = MPI_Recv(fbuf->header, bufsz, MPI_BYTE, src, HEADER_TAG, gl_comps[i].comm, &status);
//                    CHECK_MPI(errno);
//
//                    receive_data(fbuf, src, gl_comps[i].comm);
//                    fbuf->distr_ndone++;
//                    if(fbuf->distr_ndone == fbuf->distr_nranks)
//                        fbuf->is_ready = 1;
//                } else
//                    fbuf->is_ready = 1;
//
//                break;
//            case RECV_READY_TAG:
//                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].comm, &status);
//
//                CHECK_MPI(errno);
//                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//                assert(fbuf != NULL);
//                DTF_DBG(VERBOSE_DBG_LEVEL, "Receive RECV_READY_TAG notif for %s", filename);
//                Send the hearder
//                errno = MPI_Send(fbuf->header, (int)fbuf->hdr_sz, MPI_BYTE, src, HEADER_TAG, gl_comps[i].comm);
//                CHECK_MPI(errno);
//                send_data(fbuf, src, gl_comps[i].comm);
//
//                fbuf->distr_ndone++;
//                break;
//            case ROOT_MST_TAG:
//                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, ROOT_MST_TAG, gl_comps[i].comm, &status);
//                CHECK_MPI(errno);
//                DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive ROOT_MST_TAG notif for %s", filename);
//
//                fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//                assert(fbuf != NULL);
//                assert(fbuf->root_writer == -1);
//                fbuf->root_writer = src;
//                break;
//            default:
//                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: unknown tag %d", status.MPI_TAG);
//                assert(0);
//        }
//    }
//}

MPI_Offset last_1d_index(int ndims, const MPI_Offset *shape)
{
    MPI_Offset ret = 0;

    if(ndims > 0){
        int i;
        MPI_Offset *tmp = (MPI_Offset*)dtf_malloc(ndims*sizeof(MPI_Offset));
        assert(tmp != NULL);
        for(i = 0; i < ndims; i++)
            tmp[i] = shape[i] - 1;
        ret = to_1d_index(ndims, shape, tmp);
        dtf_free(tmp, ndims*sizeof(MPI_Offset));
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
    int errno;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Inside notify_file_ready");
    //comp_id = fbuf->reader_id;
    if(fbuf->iomode == DTF_IO_MODE_FILE){
        //int rank;
       // MPI_Comm_rank(fbuf->comm, &rank);

       // MPI_Comm_remote_size(gl_comps[comp_id].comm, &nranks);
        if(gl_my_rank == fbuf->root_writer){
//            if(fbuf->root_reader == -1){
//                DTF_DBG(VERBOSE_DBG_LEVEL, "Start waiting to know who is the root reader for file %s", fbuf->file_path);
//                while(fbuf->root_reader == -1)
//                    progress_io();
//            }
            if(fbuf->root_reader != -1){
                DTF_DBG(VERBOSE_DBG_LEVEL,   "Notify reader root rank %d that file %s is ready", fbuf->root_reader, fbuf->file_path);
                errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, FILE_READY_TAG, gl_comps[fbuf->reader_id].comm);
                CHECK_MPI(errno);
                fbuf->fready_notify_flag = RDR_NOTIFIED;
            }
        }
    }
}

void close_file(file_buffer_t *fbuf)
{
    /*Note: there was a barrier in the upper function*/
    DTF_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", fbuf->file_path);

    if(fbuf->writer_id == gl_my_comp_id){

        if(fbuf->iomode == DTF_IO_MODE_FILE) {
            assert(fbuf->comm != MPI_COMM_NULL);
            if(gl_my_rank == fbuf->root_writer)
                notify_file_ready(fbuf);
        } else if((fbuf->iomode == DTF_IO_MODE_MEMORY) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH)){
            /* cannot close file untill the reader closed it as well.
             then we can clean up all the requests etc. */
            DTF_DBG(VERBOSE_DBG_LEVEL, "Reader hasn't closed the file yet. Waiting...");
            while( !fbuf->rdr_closed_flag)
                progress_io_matching();
            DTF_DBG(VERBOSE_DBG_LEVEL, "Cleaning up everything");
            delete_ioreqs(fbuf);
            if(fbuf->mst_info->iodb != NULL){
                clean_iodb(fbuf->mst_info->iodb);
                dtf_free(fbuf->mst_info->iodb, sizeof(ioreq_db_t));
                fbuf->mst_info->iodb = NULL;
            }
        }
    } else if (fbuf->reader_id == gl_my_comp_id){
        if((fbuf->iomode == DTF_IO_MODE_MEMORY) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH)){
            int rank;
            assert(fbuf->rreq_cnt == 0);
            assert(fbuf->wreq_cnt == 0);
            MPI_Comm_rank(fbuf->comm, &rank);
            if(rank == 0){
                /*Notify the root writer I am closing the file*/
                errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->root_writer, IO_CLOSE_FILE_TAG, gl_comps[fbuf->writer_id].comm);
                CHECK_MPI(errno);
            }
//            /*Reader never needs these flags but set them just in case*/
//            fbuf->fclosed_flag = 1;
//            fbuf->fclose_notify_flag = 1;
        }
    }
}

void open_file(file_buffer_t *fbuf, MPI_Comm comm)
{
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Enter dtf_open %s", fbuf->file_path);

    MPI_Status status;
    int rank;

    if(fbuf->reader_id == gl_my_comp_id){

        if(fbuf->iomode == DTF_IO_MODE_FILE){

            MPI_Comm_rank(comm, &rank);
            if(fbuf->root_writer == -1){
                if(rank == 0){
                    int ncid;
                    /*First, find out who is the root master.
                      In this case, only need to copy the file name and root reader rank*/
                    void *buf = dtf_malloc(MAX_FILE_NAME+2*sizeof(int));
                    assert(buf != NULL);
                    memcpy(buf, fbuf->file_path, MAX_FILE_NAME);
                    memcpy((unsigned char*)buf+MAX_FILE_NAME, &gl_my_rank, sizeof(int));
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Asking writer who is the root of file %s", fbuf->file_path);
                    errno = MPI_Send(buf, (int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, 0, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm);
                    CHECK_MPI(errno);
                    dtf_free(buf, MAX_FILE_NAME+2*sizeof(int));

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Starting to wait for file info for %s", fbuf->file_path);
                    errno = MPI_Recv(&ncid, 1, MPI_INT, MPI_ANY_SOURCE, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
                    CHECK_MPI(errno);
                    fbuf->ncid = ncid;
                    fbuf->root_writer = status.MPI_SOURCE;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Root for file %s is %d", fbuf->file_path, fbuf->root_writer);
                }
                errno = MPI_Bcast(&fbuf->root_writer, 1, MPI_INT, 0, comm);
                CHECK_MPI(errno);
            }

            if(rank == 0){
                /*root reader rank will wait until writer finishes writing the file.
                 then it will broadcast that the file is ready to everyone else*/
                while(!fbuf->is_ready)
                    progress_io_matching();
            }

            errno = MPI_Bcast(&fbuf->is_ready, 1, MPI_INT, 0, comm);
            CHECK_MPI(errno);
            assert(fbuf->is_ready == 1);

        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
            /*First, find out who is the root master*/
            if(gl_conf.distr_mode == DISTR_MODE_REQ_MATCH){
                int nranks;
                int bufsz;
                void *buf;

                if(fbuf->root_writer != -1){
                    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: We already know the masters for this file from, e.g., a previous iteration");
                    goto fn_exit;
                }
                MPI_Comm_rank(comm, &rank);
                /*Zero rank will inquire the pnetcdf header/dtf vars info/masters info
                from writer's global zero rank and then broadcast this info to other
                readers that opened the file*/
                if(rank == 0){
                    buf = dtf_malloc(MAX_FILE_NAME+2*sizeof(int));
                    assert(buf != NULL);
                    memcpy(buf, fbuf->file_path, MAX_FILE_NAME);
                    memcpy((unsigned char*)buf+MAX_FILE_NAME, &gl_my_rank, sizeof(int));
                    MPI_Comm_size(comm, &nranks);
                    memcpy((unsigned char*)buf+MAX_FILE_NAME+sizeof(int), &nranks, sizeof(int));
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Asking writer who is the root of file %s", fbuf->file_path);
                    errno = MPI_Send(buf, (int)(MAX_FILE_NAME+sizeof(int)*2), MPI_BYTE, 0, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm);
                    CHECK_MPI(errno);
                    dtf_free(buf, MAX_FILE_NAME+2*sizeof(int));
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Starting to wait for file info for %s", fbuf->file_path);
                    errno = MPI_Probe(MPI_ANY_SOURCE, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
                    CHECK_MPI(errno);
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    fbuf->root_writer = status.MPI_SOURCE;
                    buf = dtf_malloc(bufsz);
                    assert(buf != NULL);
                    errno = MPI_Recv(buf, bufsz, MPI_BYTE, fbuf->root_writer, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
                    CHECK_MPI(errno);
                }
               // MPI_Barrier(comm);
                DTF_DBG(VERBOSE_DBG_LEVEL, "Bcast file info to others");
                errno = MPI_Bcast(&bufsz, 1, MPI_INT, 0, comm);
                CHECK_MPI(errno);
                assert(bufsz > 0);

                if(rank != 0){
                    buf = dtf_malloc(bufsz);
                    assert(buf != NULL);
                }
                errno = MPI_Bcast(buf, bufsz, MPI_BYTE, 0, comm);
                CHECK_MPI(errno);

                unpack_file_info(bufsz, buf);
                dtf_free(buf, bufsz);
                fbuf->is_ready = 1;
            }
        }
    } else if(fbuf->writer_id == gl_my_comp_id){
        //do we need it?
        assert(0);
    }
fn_exit:
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Exit dtf_open %s", fbuf->file_path);
}

int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    dtf_var_t *var = new_var(varid, ndims, dtype, shape);
    add_var(&fbuf->vars, var);
    fbuf->var_cnt++;
    return 0;
}

/*Write pnetcdf header*/
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header)
{
    DTF_DBG(VERBOSE_DBG_LEVEL, "Writing header (sz %d)", (int)hdr_sz);
    fbuf->hdr_sz = hdr_sz;
    fbuf->header = dtf_malloc(hdr_sz);
    assert(fbuf->header != NULL);
    memcpy(fbuf->header, header, (size_t)hdr_sz);
    return;
}

/*Read pnetcdf header*/
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(offset+chunk_sz > fbuf->hdr_sz){
        DTF_DBG(VERBOSE_ALL_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, fbuf->hdr_sz);
        chunk_sz = fbuf->hdr_sz - offset;
    }

    memcpy(chunk, (unsigned char*)fbuf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}


int mpitype2int(MPI_Datatype dtype)
{

    if(dtype == MPI_SIGNED_CHAR)     return 1;
    if(dtype == MPI_CHAR)            return 2;
    if(dtype == MPI_SHORT)           return 3;
    if(dtype == MPI_INT)             return 4;
    if(dtype == MPI_FLOAT)           return 5;
    if(dtype == MPI_DOUBLE)          return 6;
    if(dtype == MPI_UNSIGNED_CHAR)   return 7;
    if(dtype == MPI_UNSIGNED_SHORT)  return 8;
    if(dtype == MPI_UNSIGNED)        return 9;
    if(dtype == MPI_LONG_LONG_INT)   return 10;
    if(dtype == MPI_UNSIGNED_LONG_LONG) return 11;

    DTF_DBG(VERBOSE_ERROR_LEVEL, "Unknown mpi type");
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    return 0;
}

MPI_Datatype int2mpitype(int num)
{
    switch(num){
        case 1 :      return MPI_SIGNED_CHAR;
        case 2 :      return MPI_CHAR;
        case 3 :      return MPI_SHORT;
        case 4 :      return MPI_INT;
        case 5 :      return MPI_FLOAT;
        case 6 :      return MPI_DOUBLE;
        case 7 :      return MPI_UNSIGNED_CHAR;
        case 8 :      return MPI_UNSIGNED_SHORT;
        case 9 :      return MPI_UNSIGNED;
        case 10 :     return MPI_LONG_LONG_INT;
        case 11 :     return MPI_UNSIGNED_LONG_LONG;
        default:
            DTF_DBG(VERBOSE_ERROR_LEVEL, "Unknown mpi type");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    return MPI_DATATYPE_NULL;
}

void* dtf_malloc(size_t size)
{
    gl_stats.malloc_size += size;
    return malloc(size);
}

void dtf_free(void *ptr, size_t size)
{
    gl_stats.malloc_size -= size;
    free(ptr);
    return;
}
