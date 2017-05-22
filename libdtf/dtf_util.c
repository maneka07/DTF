/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */


#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "dtf_util.h"
#include "dtf.h"
#include "dtf_file_buffer.h"
#include "dtf_common.h"
#include "dtf_req_match.h"

extern file_info_req_q_t *gl_finfo_req_q;

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
            DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer 22, process the file info req for file %s", req->filename);
            memcpy(&fbuf->root_reader, (unsigned char*)(req->buf)+MAX_FILE_NAME, sizeof(int));

            if(fbuf->iomode == DTF_IO_MODE_MEMORY && gl_conf.distr_mode == DISTR_MODE_REQ_MATCH){
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
            errno = MPI_Send(req->buf, (int)(MAX_FILE_NAME+sizeof(int)), MPI_BYTE, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm);
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
        dtf_free(req->buf, MAX_FILE_NAME+sizeof(int));
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


MPI_Offset to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord)
{
      int i, j;
      MPI_Offset idx = 0, mem=0;

      if(ndims == 0) //scalar
        return 0;
      else if(ndims == 1) //1d array
        return *coord - *block_start;

      for(i = 0; i < ndims; i++){
        mem = coord[i] - block_start[i];
        for(j = i+1; j < ndims; j++)
          mem *= block_count[j];
        idx += mem;
      }
    return idx;
}

void notify_file_ready(file_buffer_t *fbuf)
{
    int errno;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Inside notify_file_ready");
    if(fbuf->iomode == DTF_IO_MODE_FILE){
        if(gl_my_rank == fbuf->root_writer){

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
            //reset flags
            fbuf->rdr_closed_flag = 0;

            /*File is opened and closed multiple times in SCALE-LETKF
              but it's the same set of processes, hence, don't delete the data.
            */
            /*dtf_free(fbuf->mst_info->masters, fbuf->mst_info->nmasters*sizeof(int));
            dtf_free(fbuf->mst_info, sizeof(master_info_t));
            fbuf->mst_info = NULL;
            fbuf->root_reader = -1;
            fbuf->root_writer = -1;*/
        }
    } else if (fbuf->reader_id == gl_my_comp_id){
        //if((fbuf->iomode == DTF_IO_MODE_MEMORY) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH)){

            assert(fbuf->rreq_cnt == 0);
            assert(fbuf->wreq_cnt == 0);

            if(fbuf->root_reader == gl_my_rank){
                int i;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer masters readers have closed the file");
                for(i = 0; i < fbuf->mst_info->nmasters; i++){
                    /*Notify the root writer I am closing the file*/
                    errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->mst_info->masters[i], IO_CLOSE_FILE_TAG, gl_comps[fbuf->writer_id].comm);
                    CHECK_MPI(errno);
                }
            }
            /*dtf_free(fbuf->mst_info->masters, fbuf->mst_info->nmasters*sizeof(int));
            dtf_free(fbuf->mst_info, sizeof(master_info_t));
            fbuf->mst_info = NULL;
            fbuf->root_reader = -1;
            fbuf->root_writer = -1;*/
//            /*Reader never needs these flags but set them just in case*/
//            fbuf->fclosed_flag = 1;
//            fbuf->fclose_notify_flag = 1;
        //}
    }
}

void open_file(file_buffer_t *fbuf, MPI_Comm comm)
{
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Enter dtf_open %s", fbuf->file_path);

    MPI_Status status;
    int rank;

    if(fbuf->reader_id == gl_my_comp_id){
        MPI_Comm_rank(comm, &rank);
        if(rank == 0)
            fbuf->root_reader = gl_my_rank;
        errno = MPI_Bcast(&fbuf->root_reader, 1, MPI_INT, 0, comm);
        CHECK_MPI(errno);

        if(fbuf->iomode == DTF_IO_MODE_FILE){
            if(fbuf->root_writer == -1){
                if(rank == 0){
                    int ncid;
                    MPI_Request req;
                    /*First, find out who is the root master.
                      In this case, only need to copy the file name and root reader rank*/
                    void *buf = dtf_malloc(MAX_FILE_NAME+sizeof(int));
                    assert(buf != NULL);
                    memcpy(buf, fbuf->file_path, MAX_FILE_NAME);
                    memcpy((unsigned char*)buf+MAX_FILE_NAME, &gl_my_rank, sizeof(int));
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Asking writer who is the root of file %s", fbuf->file_path);
                    errno = MPI_Isend(buf, (int)(MAX_FILE_NAME+sizeof(int)), MPI_BYTE, 0, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm, &req);
                    CHECK_MPI(errno);
                    errno = MPI_Wait(&req, MPI_STATUS_IGNORE);
                    CHECK_MPI(errno);
                    dtf_free(buf, MAX_FILE_NAME+2*sizeof(int));

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Starting to wait for file info for %s", fbuf->file_path);
                    errno = MPI_Recv(&ncid, 1, MPI_INT, MPI_ANY_SOURCE, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
                    CHECK_MPI(errno);
                    fbuf->ncid = ncid;
                    fbuf->root_writer = status.MPI_SOURCE;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Root for file %s is %d", fbuf->file_path, fbuf->root_writer);
                }
                DTF_DBG(VERBOSE_DBG_LEVEL, "Broadcast root & ncid to other readers");
                errno = MPI_Bcast(&fbuf->ncid, 1, MPI_INT, 0, comm);
                CHECK_MPI(errno);
                errno = MPI_Bcast(&fbuf->root_writer, 1, MPI_INT, 0, comm);
                CHECK_MPI(errno);
            }

            fbuf->is_ready = 0;
            if(rank == 0){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Waiting for file to become ready");
                /*root reader rank will wait until writer finishes writing the file.
                 then it will broadcast that the file is ready to everyone else*/
                double t_start = MPI_Wtime();
                while(!fbuf->is_ready)
                    progress_io_matching();
                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: Waiting to open file %.3f", MPI_Wtime()-t_start);
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

                if(fbuf->root_writer != -1)
                    goto fn_exit;   //already got all info from previous iteration

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
                    errno = MPI_Send(buf, (int)(MAX_FILE_NAME+sizeof(int)), MPI_BYTE, 0, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm);
                    CHECK_MPI(errno);
                    dtf_free(buf, MAX_FILE_NAME+sizeof(int));
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

        //Notify writer
        if(rank == 0){
            errno = MPI_Send(&fbuf->ncid, 1, MPI_INT, fbuf->root_writer, IO_OPEN_FILE_FLAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(errno);
        }
    } else if(fbuf->writer_id == gl_my_comp_id){
        //do we need it?
        assert(0);
        /*reset all flags*/
    }
fn_exit:
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Exit dtf_open %s", fbuf->file_path);
}

int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int el_sz;
    int i;
    dtf_var_t *var = new_var(varid, ndims, dtype, shape);
    add_var(fbuf, var);
    MPI_Type_size(dtype, &el_sz);

    DTF_DBG(VERBOSE_DBG_LEVEL, "varid %d, dim %d, el_sz %d. shape:", varid, ndims, el_sz);
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t%lld", shape[i]);

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

/*Find the biggest subblock of data that can fit into given buffer size sbufsz
  [INOUT] cur_count - the subblock that fits
  [INOUT] cur_nelems - how many elements managed to fit. Zero means didn't fit anything.
  [IN] the rest
*/
void find_fit_block(int ndims,
		    int cur_dim,
		    const MPI_Offset *start,
		    const MPI_Offset *count,
		    MPI_Offset *cur_start,
		    MPI_Offset *cur_count,
		    const size_t sbufsz,
		    const size_t el_sz,
		    MPI_Offset *cur_nelems,
		    MPI_Offset tot_nelems)
{
  int i;

  if(*cur_nelems == tot_nelems)
    return;
  if(sbufsz == 0)
    return;

  /*Compute if current subblock can fit*/
  MPI_Offset full_subbl = 1;
  for(i = 0; i < ndims; i++)
    if(i != cur_dim)
      full_subbl *= cur_count[i];

  DTF_DBG(VERBOSE_DBG_LEVEL, "full subblock %lld, cur dim %d, ndims %d", full_subbl, cur_dim, ndims);

  for(i = count[cur_dim] - (cur_start[cur_dim] - start[cur_dim]); i > 0; i--){
    //DTF_DBG(VERBOSE_DBG_LEVEL, "left %lld, right %lld", full_subbl * i * el_sz, (MPI_Offset)sbufsz);
    if(full_subbl * i * el_sz <= (MPI_Offset)sbufsz)
      break;
  }

  cur_count[cur_dim] = i;
  assert(cur_count[cur_dim] > 0);

  *cur_nelems = 1;
  for(i = 0; i < ndims; i++){
    *cur_nelems *= cur_count[i];
  }

  if(cur_count[cur_dim] == count[cur_dim]){
    if(cur_dim == 0){
      assert(*cur_nelems == tot_nelems);
      return;
    } else {
      DTF_DBG(VERBOSE_DBG_LEVEL, "Go higher. Cur nelems %lld. Cur count:", *cur_nelems);
      for(i=0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "  %lld", cur_count[i]);
      /*Go higher one dimension*/
      find_fit_block(ndims, cur_dim - 1, start, count, cur_start, cur_count, sbufsz, el_sz, cur_nelems, tot_nelems);
    }
  }
}

void shift_coord(int ndims, const MPI_Offset *bl_start,
                 const MPI_Offset *bl_count, MPI_Offset *subbl_start,
                 MPI_Offset *subbl_count, MPI_Offset fit_nelems)
{
    int i;

    /*Shift the start position*/
    if(fit_nelems == 1){ //special case
      subbl_start[ndims-1]++;
    } else {
        for(i = 0; i < ndims; i++)
            if(subbl_count[i] > 1)
                  subbl_start[i] += subbl_count[i];
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "New start before adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t %lld", subbl_start[i]);

    for(i = ndims - 1; i > 0; i--)
        if(subbl_start[i] == bl_start[i] + bl_count[i]){
            subbl_start[i] = bl_start[i];
            if( (subbl_start[i-1] != bl_start[i-1] + bl_count[i-1]) && (subbl_count[i-1] == 1)){
                subbl_start[i-1]++;
            }
        } else
            break;

    DTF_DBG(VERBOSE_DBG_LEVEL, "New start after adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t %lld", subbl_start[i]);


//    DTF_DBG(VERBOSE_DBG_LEVEL, "Copied subblock. Shift start:");
//    for(i = 0; i < var->ndims; i++)
//        DTF_DBG(VERBOSE_DBG_LEVEL, "   %lld\t -->\t %lld", bl_start[i], subbl_start[i]);
}

double compute_checksum(void *arr, int ndims, const MPI_Offset *shape, MPI_Datatype dtype)
{
    double sum = 0;
    unsigned nelems;
    int i;

    if(dtype != MPI_DOUBLE && dtype != MPI_FLOAT){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: checksum supported only for double or float data");
        return 0;
    }

    if(ndims == 0){
        if(dtype == MPI_DOUBLE)
            sum = *(double*)arr;
        else
            sum = *(float*)arr;
        return sum;
    }

    nelems = shape[0];
    for(i = 1; i < ndims; i++)
        nelems *= shape[i];

    for(i = 0; i < nelems; i++)
        if(dtype == MPI_DOUBLE)
            sum += ((double*)arr)[i];
        else
            sum += ((float*)arr)[i];
    return sum;
}

/*only support conversion double<->float*/
void get_put_data(dtf_var_t *var,
                  MPI_Datatype dtype,
                  unsigned char *block_data,
                  const MPI_Offset *block_start,
                  const MPI_Offset *block_count,
                  const MPI_Offset subbl_start[],
                  const MPI_Offset subbl_count[],
                  unsigned char *subbl_data,
                  int get_put_flag,
                  int convert_flag)
{
    int i;
    MPI_Offset *cur_coord = dtf_malloc(var->ndims*sizeof(MPI_Offset));
    for(i = 0; i < var->ndims; i++){
        cur_coord[i] = subbl_start[i];
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Call getput data");
    /*read from user buffer to send buffer*/
    recur_get_put_data(var, dtype, block_data, block_start,
                       block_count, subbl_start, subbl_count, 0,
                       cur_coord, subbl_data, get_put_flag, convert_flag);
    dtf_free(cur_coord, var->ndims*sizeof(MPI_Offset));
    gl_stats.ngetputcall++;
}

/*only support conversion double<->float*/
void recur_get_put_data(dtf_var_t *var,
                          MPI_Datatype dtype,
                          unsigned char *block_data,
                          const MPI_Offset *block_start,
                          const MPI_Offset *block_count,
                          const MPI_Offset subbl_start[],
                          const MPI_Offset subbl_count[],
                          int dim,
                          MPI_Offset coord[],
                          unsigned char *subbl_data,
                          int get_put_flag,
                          int convert_flag)
{
    int i;
    if(dim == var->ndims - 1){
        int bl_type_sz, subbl_type_sz;
        MPI_Type_size(dtype, &bl_type_sz);
        MPI_Type_size(var->dtype, &subbl_type_sz);
        MPI_Offset block_offt = to_1d_index(var->ndims, block_start, block_count, coord)*bl_type_sz;
        MPI_Offset subbl_offt = to_1d_index(var->ndims, subbl_start, subbl_count, coord)*subbl_type_sz;
        MPI_Offset nelems = subbl_count[var->ndims-1];
        //MPI_Offset data_sz = subbl_count[var->ndims-1]*subbl_type_sz;    //data_sz

        if(get_put_flag == DTF_READ){
            if(convert_flag){
                convertcpy(dtype, var->dtype, (void*)(block_data+block_offt), (void*)(subbl_data+subbl_offt), (int)nelems);
//                for(i = 0; i < subbl_count[var->ndims-1]; i++)
//                    ((float*)(subbl_data+subbl_offt))[i] = (float)((double*)(block_data+block_offt))[i];
            } else
                /*copy data block -> subblock*/
                memcpy(subbl_data+subbl_offt, block_data+block_offt, nelems*subbl_type_sz);
        } else { /*DTF_WRITE*/
           /*copy data subblock -> block*/
            if(convert_flag)
                convertcpy(var->dtype, dtype, (void*)(subbl_data+subbl_offt),(void*)(block_data+block_offt), nelems);
            else
                memcpy(block_data+block_offt, subbl_data+subbl_offt, nelems*subbl_type_sz);
        }
        return;
    }

    for(i = 0; i < subbl_count[dim]; i++){
        coord[dim] = subbl_start[dim] + i;
        recur_get_put_data(var, dtype, block_data, block_start, block_count,
                           subbl_start, subbl_count, dim+1, coord, subbl_data,
                           get_put_flag, convert_flag);
    }
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
    if(size > gl_stats.malloc_size)
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: mem stat negative (left %lu), will free %lu", gl_stats.malloc_size, size);
    gl_stats.malloc_size -= size;
    free(ptr);
    return;
}

void convertcpy(MPI_Datatype type1, MPI_Datatype type2, void* srcbuf, void* dstbuf, int nelems)
{
    int i;
    if(type1 == MPI_FLOAT){
        assert(type2 == MPI_DOUBLE);
        for(i = 0; i < nelems; i++)
            ((double*)dstbuf)[i] = (double)(((float*)srcbuf)[i]);
    } else if(type1 == MPI_DOUBLE){
        assert(type2 = MPI_FLOAT);
        for(i = 0; i < nelems; i++)
            ((float*)dstbuf)[i] = (float)(((double*)srcbuf)[i]);
    }
}
