/*
*
* Buffered I/O:  we internally buffer data on the writer's side and on the reader's side
*
*
*/
#include <assert.h>
#include "pfarb_buf_io.h"
#include "pfarb_util.h"
#include "pfarb.h"
#include "pfarb_mem.h"
#include "pfarb_common.h"

//TODO remove addressing file by name to addressing by ncid
MPI_Offset buf_read_write_var( file_buffer_t *fbuf,
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
    int el_sz;
    if(rw_flag == FARB_READ)
         assert(fbuf->is_ready);

    if(imap != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: could not find var with id %d", varid);
        return 0;
    }

    /*check number of elements to read*/
    if(count != NULL){
        MPI_Offset nelems;
        int i;
        nelems = count[0];
        for(i = 1; i < var->ndims; i++)
            nelems *= count[i];

        if(nelems == 0){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Nothing to read or write");
            return 0;
        }

    }

    MPI_Type_size(dtype, &el_sz);
    assert(el_sz > 0);
    FARB_DBG(VERBOSE_ERROR_LEVEL, "Mismatch of var's datatype");
    assert(var->dtype == dtype);
    //TODO implement type conversion

    if(var->ndims <= 1){ /*scalar or 1d array*/
        MPI_Offset offt = (*start)*el_sz;
        MPI_Offset data_sz = (*count)*el_sz;

        if(rw_flag == FARB_READ){
            ret = mem_contiguous_read(var, offt, data_sz, buf);
        } else {
            ret = mem_contiguous_write(var,offt, data_sz, buf);
        }

        if(ret != data_sz)
            FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: Meant to read/write %llu bytes but actual val is %llu", data_sz, ret);

    } else { /*multi-dimentional array*/
        if(rw_flag == FARB_WRITE){
            ret = mem_noncontig_write(var, start, count, buf);
        } else { /*FARB_READ*/
            ret = mem_noncontig_read(var, start, count, buf);
        }
    }

    return ret;
}


/*
    vars        [OUT]
    var_cnt     [OUT]
    the rest    [IN]
*/
static void unpack_vars(file_buffer_t *fbuf, int buf_sz, void *buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord)
{
    size_t offt = 0;
    MPI_Offset offset, data_sz, ncnt;
    int var_cnt;
    assert(buf_sz > 0);
    assert(buf != NULL);
    *node_cnt = 0;
    *first_el_coord = NULL;
    int first_el_coord_sz = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    int type;

    FARB_DBG(VERBOSE_DBG_LEVEL, "Unpacking vars sz %d", buf_sz);

    int i, j, varid;
    farb_var_t *var;
    buffer_node_t *node;

    var_cnt = *((int*)chbuf);
    offt += sizeof(int);
    FARB_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", var_cnt);
    for(i = 0; i < var_cnt; i++){

        varid = *((int*)(chbuf+offt));

        offt += sizeof(int);
        var = find_var(fbuf->vars, varid);
        if(var == NULL){
            var = new_var(varid, 0, 0, NULL);
            type = (int)(*((MPI_Offset*)(chbuf+offt)));
            var->dtype = int2mpitype(type);
            offt+=sizeof(MPI_Offset);
            var->ndims = *((int*)(chbuf+offt));
            offt += sizeof(int);

            if(var->ndims > 0){
                var->shape = malloc(var->ndims*sizeof(MPI_Offset));
                assert(var->shape != NULL);
                memcpy((void*)var->shape, chbuf+offt, sizeof(MPI_Offset)*var->ndims);
                offt += sizeof(MPI_Offset)*var->ndims;
            } else
                var->shape = NULL;

            add_var(&(fbuf->vars), var);
            fbuf->var_cnt++;
        } else {
            /*Skip the fields that have already been defined*/
            offt += sizeof(MPI_Offset) + sizeof(int) + sizeof(MPI_Offset)*var->ndims;
        }

        ncnt = *((MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        /*Unpack data about nodes to receive*/
        if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
            var->node_cnt += ncnt;
            *node_cnt += ncnt;

            /*Unpack node offsets*/
            for(j = 0; j < (int)ncnt; j++){
                offset = *((MPI_Offset*)(chbuf+offt));
                offt += sizeof(MPI_Offset);
                data_sz = *((MPI_Offset*)(chbuf+offt));
                offt += sizeof(MPI_Offset);
                node = new_buffer_node(offset, data_sz, 0);
                insert_buffer_node(&var->nodes, node);
            }
        } else { /*DIST_PATTERN_SCATTER*/
            if(ncnt != 0) {
                (*node_cnt)++;
                /*Only unpack the info about distr_count[] and first_coord[]*/
                MPI_Offset *tmp = (MPI_Offset*)realloc(var->distr_count, var->ndims*sizeof(MPI_Offset));
                assert(tmp != NULL);
                var->distr_count = tmp;
                memcpy((void*)var->distr_count, chbuf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);

                first_el_coord_sz += var->ndims;
                tmp = (MPI_Offset*)realloc(*first_el_coord, first_el_coord_sz*sizeof(MPI_Offset));
                assert(tmp != NULL);
                *first_el_coord = tmp;

                // need to remember the coordinate of the first (corner) element to be able to unpack
                // contiguous 1D buffer into multi-dimensional nodes
                for(j = 0; j < var->ndims; j++){
                    (*first_el_coord)[first_el_coord_sz - var->ndims + j] = *((MPI_Offset*)(chbuf+offt));
                    offt += sizeof(MPI_Offset);
                }
            } else {
                if(var->distr_count != NULL) //could be allocated during recv from another writer rank
                    free(var->distr_count);
                var->distr_count = NULL;
            }
        }
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Unpacking vars: offt is %lu, bufsz is %lu", offt, (size_t)buf_sz);
    assert(offt == (size_t)buf_sz);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished unpacking vars");
}

int receive_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm)
{
    MPI_Status status;
    int buf_sz;
    int errno, i;
    void *buf;
    farb_var_t *var;
    buffer_node_t *node;
    MPI_Offset node_cnt;
    int ncnt = 0;
    MPI_Offset *first_el_coord;

    /*Receive vars*/
    MPI_Probe(rank, VARS_TAG, intercomm, &status);
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &buf_sz);
    buf = malloc((size_t)buf_sz);
    assert(buf != NULL);
    errno = MPI_Recv(buf, buf_sz, MPI_UNSIGNED_CHAR, rank, VARS_TAG, intercomm, &status);
    CHECK_MPI(errno);

    unpack_vars(fbuf, buf_sz, buf, &node_cnt, &first_el_coord);
    free(buf);
    /*Receive memory nodes*/
    var = fbuf->vars;

    if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will recv %d nodes", (int)node_cnt);
        MPI_Request *reqs = NULL;
        reqs = (MPI_Request*)malloc(sizeof(MPI_Request)*node_cnt);
        assert(reqs != NULL);
        ncnt = 0;

        while(var != NULL){
            print_nodes(var->nodes);

            node = var->nodes;
            while(node != NULL){
                if(node->data != NULL){
                    node = node->next;
                    continue;
                }
                FARB_DBG(VERBOSE_ALL_LEVEL, "Receive node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                node->data = malloc((size_t)node->data_sz);
                assert(node->data != NULL);
                errno = MPI_Irecv(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm, &reqs[ncnt]);
                CHECK_MPI(errno);
                ncnt++;
                node = node->next;
                if(ncnt == (int)node_cnt)
                    break;
            }
            var = var->next;
            if(ncnt == (int)node_cnt)
                break;
        }
        errno = MPI_Waitall(node_cnt, reqs, MPI_STATUSES_IGNORE);
        CHECK_MPI(errno);
        free(reqs);
    } else { /*DISTR_PATTERN_SCATTER*/
        void *rbuf = NULL, *tmp;
        int bufsz;
        int first_el_coord_sz = 0;
        MPI_Offset written;
        int el_sz;

        while(var != NULL){
            if(var->distr_count == NULL){
                var = var->next;
                continue;
            }
            bufsz = 1;
            for(i = 0; i < var->ndims; i++)
                bufsz *= (int)var->distr_count[i];
            MPI_Type_size(var->dtype, &el_sz);
            bufsz *= el_sz;
            tmp = realloc(rbuf, bufsz);
            assert(tmp != NULL);
            rbuf = tmp;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Try to recv node of sz %d", (int)bufsz);
            errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm, MPI_STATUS_IGNORE);
            CHECK_MPI(errno);

            written = mem_noncontig_write(var, &first_el_coord[first_el_coord_sz], var->distr_count, rbuf);
            assert((int)written == bufsz);
            first_el_coord_sz += (int)var->ndims;
            ncnt++;
            var = var->next;
        }
        free(rbuf);
    }
    assert(ncnt == (int)node_cnt);
    return 0;
}


static void pack_vars(file_buffer_t *fbuf, int dst_rank, int *buf_sz, void **buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord)
{
    int i;
    size_t sz = 0, offt=0;
    buffer_node_t *node;
    sz += sizeof(fbuf->var_cnt);
    *node_cnt = 0;
    *first_el_coord = NULL;
    int first_el_coord_sz = 0;
    unsigned char* chbuf;
    farb_var_t *var = fbuf->vars;
    while(var != NULL){
        sz += sizeof(var->id) + sizeof(MPI_Offset)/*dtype*/ + sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims +
        sizeof(MPI_Offset) + var->node_cnt*sizeof(MPI_Offset)*2;
        var = var->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Packing vars: sz %lu", sz);

    *buf = malloc(sz);
    assert(*buf != NULL);
    *node_cnt = 0;
    chbuf = (unsigned char*)(*buf);
    //write the number of vars
    *((int*)(chbuf)) = fbuf->var_cnt;
    offt += sizeof(int);

    var = fbuf->vars;
    while(var != NULL){
        *((int*)(chbuf+offt)) = var->id;
        offt += sizeof(int);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)mpitype2int(var->dtype);
        offt += sizeof(MPI_Offset);
        *((int*)(chbuf+offt)) = var->ndims;
        offt += sizeof(int);

        memcpy(chbuf+offt, (void*)var->shape, sizeof(MPI_Offset)*var->ndims);
        offt += sizeof(MPI_Offset)*var->ndims;

        if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
            *((MPI_Offset*)(chbuf+offt)) = var->node_cnt;
            offt+=sizeof(MPI_Offset);
            *node_cnt += var->node_cnt;
            /*Pack node offsets*/
            node = var->nodes;
            while(node != NULL){
                *((MPI_Offset*)(chbuf+offt)) = node->offset;
                offt+=sizeof(MPI_Offset);
                *((MPI_Offset*)(chbuf+offt)) = node->data_sz;
                offt+=sizeof(MPI_Offset);
                node = node->next;
            }
        } else { /*DIST_PATTERN_SCATTER*/
            /* We will distribute blocks of count[] elemenets one by one
               to each process in the range. If we fail to read
               exactly count[] elements, it means some values are
               lacking. This will cause the app to abort.
               Writer process may have
               written the blocks non-contiguously, so we need
               to find out the start element of the block that will be sent
               to given dst_rank.
               //TODO figure out with fill values in pnetcdf
            */

            //TODO what if it's a scalar var and we need to send it?
            if(var->distr_count == NULL){
                /*If distr_count wasn't set we won't be distributed this variable to
                  the reader. But print out warning just in case if user forgot
                  to set distr_count for this variable */
                  FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: farb_set_distr_count was not called for variable with id %d. \
                  This variable will not be distributed.", var->id);
                  *((MPI_Offset*)(chbuf+offt)) = 0;
                  offt+=sizeof(MPI_Offset);
            } else {
                int j;
                *((MPI_Offset*)(chbuf+offt)) = 1;
                offt+=sizeof(MPI_Offset);
                (*node_cnt)++;
                int el_sz;

                MPI_Offset offset, index_1d;
                int rank_idx;
                MPI_Offset last_el_idx = last_1d_index(var->ndims, var->shape);

                MPI_Type_size(var->dtype, &el_sz);

                first_el_coord_sz += var->ndims;
                MPI_Offset *coord = (MPI_Offset*)realloc(*first_el_coord, first_el_coord_sz*sizeof(MPI_Offset));
                assert(coord != NULL);
                *first_el_coord = coord;

                coord = &( (*first_el_coord)[first_el_coord_sz - var->ndims] );
                /*First, find out the coordinate of the start element to send*/
                memcpy((void*)coord, (void*)var->first_coord, var->ndims*sizeof(MPI_Offset));

                for(j = 0; j < var->ndims; j++){
                    FARB_DBG(VERBOSE_ALL_LEVEL, "start with coord %d", (int)coord[j]);
                }
                offset = var->nodes->offset;
                for(rank_idx = 0; rank_idx < fbuf->distr_nranks; rank_idx++){
                    while(1){
                        index_1d = to_1d_index(var->ndims, var->shape, coord);
                        FARB_DBG(VERBOSE_ALL_LEVEL, "1d idx %d", (int)index_1d);
                        /*Check that we haven' t gone out of bounds*/
                        if(index_1d > last_el_idx){
                            FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: There is no data to send to rank %d. Incorrect config file?", dst_rank);
                            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
                        }
                        offset = (MPI_Offset)el_sz * index_1d;
                        FARB_DBG(VERBOSE_ALL_LEVEL, "offset %d", (int)offset);
                        if(data_present_at_offt(var->nodes, offset))
                            break;

                        /*find next offset at which data is present*/
                        for(i = 0; i < var->ndims; i++)
                            coord[i] += var->distr_count[i];
                    }

                    if(fbuf->distr_ranks[rank_idx] == dst_rank)
                        break;
                }

                for(j = 0; j < var->ndims; j++){
                    FARB_DBG(VERBOSE_ALL_LEVEL, "data starts at %d, count %d, shape %d", (int)coord[j], (int)var->distr_count[j], (int)var->shape[j]);
                }
                FARB_DBG(VERBOSE_DBG_LEVEL, "Will send to rank %d data starting at idx %d", dst_rank, (int)index_1d);
                assert(rank_idx != fbuf->distr_nranks);
                /*Will eventually pack multi-dimensional data to one contiguous
                  buffer to send.
                  Save the distr_count[] and first_coord[] values so that
                  the reader could unpack the data correctly later*/

                for(i=0; i < var->ndims; i++){
                    *((MPI_Offset*)(chbuf+offt)) = var->distr_count[i];
                    offt+=sizeof(MPI_Offset);
                }
                for(i=0; i < var->ndims; i++){
                    *((MPI_Offset*)(chbuf+offt)) = coord[i];
                    offt+=sizeof(MPI_Offset);
                }
            }
        }
        var = var->next;
    }

    //assert(offt == sz);

    *buf_sz = (int)offt;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Actually packed %d", *buf_sz);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finish packing vars");
}


/*Set how many elements in each dimension to distribute to ranks.
  File's distrib_pattern must be set to scatter */
int set_distr_count(file_buffer_t *fbuf, int varid, int count[])
{
    int i;
    MPI_Offset *cnt;
    if(fbuf->distr_pattern != DISTR_PATTERN_SCATTER){
        FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: cannot set distribute count. File's distribute pattern must be <scatter>.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: set_distr_count(): No variable with such id (%d). Ignored.", varid);
        return 1;
    }

    cnt = (MPI_Offset*)realloc(var->distr_count, var->ndims*sizeof(MPI_Offset));
    assert(cnt != NULL);
    var->distr_count = cnt;
    for(i = 0; i < var->ndims; i++)
        var->distr_count[i] = (MPI_Offset)count[i];

    return 0;
}


int send_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm)
{
    int errno,i;
    int buf_sz;
    void *buf;
    farb_var_t *var;
    buffer_node_t *node;
    MPI_Offset node_count;
    int ncnt = 0;
    int first_el_coord_sz = 0;
    MPI_Offset *first_el_coord;

    /*Pack the vars info and send it*/
    pack_vars(fbuf, rank, &buf_sz, &buf, &node_count, &first_el_coord);
    errno = MPI_Send(buf, buf_sz, MPI_BYTE, rank, VARS_TAG, intercomm);
    CHECK_MPI(errno);

    /*Send buffer nodes*/
    var = fbuf->vars;

    if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
        MPI_Request *reqs = NULL, *tmp;
        while(var != NULL){
            tmp = (MPI_Request*)realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
            assert(tmp != NULL);
            reqs = tmp;

            node = var->nodes;
            FARB_DBG(VERBOSE_DBG_LEVEL, "Will send %d nodes", (int)var->node_cnt);
            for(i = 0; i < var->node_cnt; i++){
                FARB_DBG(VERBOSE_ALL_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                errno = MPI_Isend(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG + ncnt, intercomm, &reqs[i]);
                CHECK_MPI(errno);
                node = node->next;
                ncnt++;
            }
            errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
            CHECK_MPI(errno);
            var = var->next;
        }
        free(reqs);
    } else { /*SCATTER*/
        void *sbuf = NULL, *tmp;
        int bufsz;
        MPI_Offset readsz;
        int el_sz;

        while(var != NULL){
            if(var->distr_count == NULL){
                var = var->next;
                continue;
            }
            //pack multi-dim data into 1d contiguous buffer
            MPI_Type_size(var->dtype, &el_sz);
            bufsz = 1;
            for(i = 0; i < var->ndims; i++)
                bufsz *= (int)var->distr_count[i];
            bufsz *= el_sz;

            tmp = realloc(sbuf, bufsz);
            assert(tmp != NULL);
            sbuf = tmp;
            for(i = 0; i < var->ndims; i++){
                FARB_DBG(VERBOSE_ALL_LEVEL, "send data starts at %d", (int)first_el_coord[first_el_coord_sz+i]);
            }
            readsz = mem_noncontig_read(var, &first_el_coord[first_el_coord_sz], var->distr_count, sbuf);
            if( (int)readsz != bufsz ){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: not enough data to send to rank %d", rank);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Should read %d but read %d", (int)bufsz, (int)readsz);
            }
            //if(readsz != 0){
            FARB_DBG(VERBOSE_ALL_LEVEL, "Sending node of sz %d", (int)readsz);
            errno = MPI_Send(sbuf, (int)readsz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm);
            CHECK_MPI(errno);
            ncnt++;
            first_el_coord_sz += var->ndims;
            //}
            var = var->next;
        }

        free(sbuf);
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "Total sent mem nodes %d, should have sent %d", ncnt, (int)node_count);
    assert((int)node_count == ncnt);
    return 0;
}




