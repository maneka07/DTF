#include "dtf_req_match.h"
#include "dtf_util.h"
#include "dtf.h"
#include <unistd.h>

file_info_req_q_t *gl_finfo_req_q = NULL;
dtf_msg_t *gl_msg_q = NULL;

static void send_data_wrt2rdr(void* buf, int bufsz);

void init_iodb(file_buffer_t *fbuf)
{
    fbuf->mst_info->iodb = dtf_malloc(sizeof(struct ioreq_db));
    assert(fbuf->mst_info->iodb != NULL);
    fbuf->mst_info->iodb->ritems = NULL;
    fbuf->mst_info->iodb->witems = NULL;
    fbuf->mst_info->iodb->nritems = 0;
    fbuf->mst_info->iodb->updated_flag = 0;
}

io_req_t *find_io_req(io_req_t *list, int var_id)
{
    io_req_t *ioreq = list;
    while(ioreq != NULL){
        if(ioreq->var_id == var_id)
            break;
        ioreq = ioreq->next;
    }
    return ioreq;
}

void delete_ioreq(file_buffer_t *fbuf, io_req_t **ioreq);
void delete_ioreqs(file_buffer_t *fbuf)
{
    io_req_t *ioreq;
    DTF_DBG(VERBOSE_ALL_LEVEL, "Delete io requests for file %s", fbuf->file_path);

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        delete_ioreq(fbuf, &ioreq);
        ioreq = fbuf->ioreqs;
    }
    assert(fbuf->rreq_cnt == 0);
    assert(fbuf->wreq_cnt == 0);
}

static void pack_file_info(file_buffer_t *fbuf, MPI_Offset *bufsz, void **buf)
{
    dtf_var_t *var;
    MPI_Offset sz = 0, offt = 0;
    unsigned char *chbuf;
    int i;
    rb_red_blk_node *var_node;
    /*Pack:
       - file name
       - file header size
       - header
       - number of masters
       - master list
       - number of vars
       - vars
    */

    sz =   MAX_FILE_NAME + fbuf->hdr_sz +
           fbuf->mst_info->nmasters*sizeof(MPI_Offset) +
           sizeof(MPI_Offset)*3;

    sz += sz%sizeof(MPI_Offset); //padding

    dtf_var_t tmp_var;
    for(i = 0; i < fbuf->var_cnt; i++){
        tmp_var.id = i;
        var_node = RBExactQuery(fbuf->vars, &tmp_var);
        assert(var_node != NULL);
        var = (dtf_var_t*)(var_node->key);
        /*sz += sizeof(var->id) + sizeof(var->el_sz) +
        sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims;*/
        sz += sizeof(MPI_Offset)*3+ sizeof(MPI_Offset)*var->ndims;
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Packing info: sz %lld", sz);
    *buf = dtf_malloc(sz);
    assert(*buf != NULL);
    chbuf = (unsigned char*)(*buf);

    /*filename*/
    memcpy(chbuf, fbuf->file_path, MAX_FILE_NAME);
    offt += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    /*header size*/
    *((MPI_Offset*)(chbuf+offt)) = fbuf->hdr_sz;
    offt += sizeof(MPI_Offset);
    /*header*/
    memcpy(chbuf+offt, fbuf->header, fbuf->hdr_sz);
    offt += fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset);
    /*number of masters*/
    *((MPI_Offset*)(chbuf+offt)) = fbuf->mst_info->nmasters;
    offt += sizeof(MPI_Offset);
    if(fbuf->mst_info->nmasters){
        DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d masters", fbuf->mst_info->nmasters);
        /*list of masters*/
        memcpy(chbuf+offt, fbuf->mst_info->masters, fbuf->mst_info->nmasters*sizeof(int));
        offt += fbuf->mst_info->nmasters*sizeof(MPI_Offset); //sizeof(int) + padding for MPI_Offset
    }

    /*number of vars*/
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->var_cnt;
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d vars", fbuf->var_cnt);
    /*vars*/

    for(i = 0; i < fbuf->var_cnt; i++){
        tmp_var.id = i;
        var_node = RBExactQuery(fbuf->vars, &tmp_var);
        assert(var_node != NULL);
        var = (dtf_var_t*)(var_node->key);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->id;
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)mpitype2int(var->dtype);
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->ndims;
        offt += sizeof(MPI_Offset);
        memcpy((void*)(chbuf+offt), (void*)var->shape, sizeof(MPI_Offset)*var->ndims);
        offt += sizeof(MPI_Offset)*var->ndims;
    }
    DTF_DBG(VERBOSE_ALL_LEVEL, "offt %lld", offt);
    assert(offt == sz);
    *bufsz = sz;
}

void unpack_file_info(MPI_Offset bufsz, void *buf)
{
    int i, varid, var_cnt;
    file_buffer_t *fbuf;
    dtf_var_t *var;
    MPI_Offset offt = 0;
    int type;
    unsigned char *chbuf = (unsigned char*)buf;
    char filename[MAX_FILE_NAME];
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start unpackinf file info for of sz %d", (int)bufsz);
    /*Unpack:
       - file name
       - file header size
       - header
       - number of masters
       - master list
       - number of vars
       - vars
    */

    /*filename*/
    memcpy(filename, chbuf, MAX_FILE_NAME);
    offt += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "unpack filename %s", filename);
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);
    /*ncid*/
    fbuf->ncid = -1;
    /*header size*/
    fbuf->hdr_sz = *((MPI_Offset*)(chbuf+offt));
    offt += sizeof(MPI_Offset);
    /*header*/
    fbuf->header = dtf_malloc(fbuf->hdr_sz);
    assert(fbuf->header != NULL);
    memcpy(fbuf->header, chbuf+offt, fbuf->hdr_sz);
    offt += fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset);
    /*number of masters*/
    fbuf->mst_info->nmasters = (int)(*((MPI_Offset*)(chbuf+offt)));
    DTF_DBG(VERBOSE_DBG_LEVEL, "unpack %d masters", fbuf->mst_info->nmasters);
    offt += sizeof(MPI_Offset);
    if(fbuf->mst_info->nmasters){
        /*list of masters*/
        fbuf->mst_info->masters = dtf_malloc(fbuf->mst_info->nmasters*sizeof(int));
        assert(fbuf->mst_info->masters != NULL);
        memcpy(fbuf->mst_info->masters, chbuf+offt, fbuf->mst_info->nmasters*sizeof(int));
        offt += fbuf->mst_info->nmasters*sizeof(MPI_Offset);
        //init root master
        fbuf->root_writer = fbuf->mst_info->masters[0];
    } else {
        fbuf->mst_info->masters = NULL;
        fbuf->root_writer = -1;
    }

    /*number of vars*/
    var_cnt = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", var_cnt);
    /*vars*/
    for(i = 0; i < var_cnt; i++){
        varid = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        var = find_var(fbuf, varid);
        assert(var == NULL);
        var = new_var(varid, 0, 0, NULL);
        type = (int)(*((MPI_Offset*)(chbuf+offt)));
        var->dtype = int2mpitype(type);
        offt+=sizeof(MPI_Offset);
        var->ndims = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);

        if(var->ndims > 0){
            var->shape = dtf_malloc(var->ndims*sizeof(MPI_Offset));
            assert(var->shape != NULL);
            memcpy((void*)var->shape, chbuf+offt, sizeof(MPI_Offset)*var->ndims);
            offt += sizeof(MPI_Offset)*var->ndims;
        } else
            var->shape = NULL;

        add_var(fbuf, var);
        fbuf->var_cnt++;
    }
    assert(offt == bufsz);
    DTF_DBG(VERBOSE_ALL_LEVEL, "Finished unpacking");
}

/*match subblocks of data*/
static void do_matching(file_buffer_t *fbuf, int intracomp_io_flag)
{
    int mlc_ranks = 4;
    int mlc_buf   = 1024;
    int i;
    write_db_item_t *witem = NULL;
    read_db_item_t *ritem;
    read_dblock_t *rblock;
    write_dblock_t *wblock;
    int rank_idx, var_id;
    int nwriters, allocd_nwriters;
    int *writers;
    unsigned char **sbuf;
    int *bufsz;
    size_t *offt;
    int matched_rank;
    MPI_Offset *matched_count;
    dtf_var_t *var = NULL;
    int ndims;
    int ntimes_while = 0;

    double t_st;

    int n_matched_blocks = 0;
    if(!fbuf->mst_info->iodb->updated_flag) //no new info since last time matching was done, ignore
        return;

    fbuf->mst_info->iodb->updated_flag = 0; //reset

    writers = (int*)dtf_malloc(mlc_ranks*sizeof(int));
    assert(writers != NULL);
    sbuf = (unsigned char**)dtf_malloc(mlc_ranks*sizeof(unsigned char*));
    assert(sbuf != NULL);
    bufsz = (int*)dtf_malloc(mlc_ranks*sizeof(int));
    assert(bufsz != NULL);
    offt = (size_t*)dtf_malloc(mlc_ranks*sizeof(size_t));
    assert(offt != NULL);

    allocd_nwriters = mlc_ranks;

    for(i = 0; i < mlc_ranks; i++){
        sbuf[i] = NULL;
        bufsz[i] = 0;
        offt[i] = 0;
    }

//    DTF_DBG(VERBOSE_ERROR_LEVEL, "before matching: %d ritems", (int)fbuf->mst_info->iodb->nritems);
    /*Try to match as many read and write
    requests as we can*/
    nwriters = 0;
    ritem = fbuf->mst_info->iodb->ritems;
    while(ritem != NULL){

        /*If we are supposed to match intra component io requests
        then skip all other requests*/
        if(intracomp_io_flag && (ritem->comm != gl_comps[gl_my_comp_id].comm)){
            ritem = ritem->next;
            continue;
        }

        if(!intracomp_io_flag && (ritem->comm == gl_comps[gl_my_comp_id].comm)){
            //for now put assert here, need to think if such situation
            //should ever happen
            DTF_DBG(VERBOSE_DBG_LEVEL, "Writer has unfinished intracomp rreqs for rank %d", ritem->rank);
            assert(0);
            ritem = ritem->next;
            continue;
        }

        n_matched_blocks = 0;

        t_st = MPI_Wtime();

        for(i = 0; i < allocd_nwriters; i++){
            sbuf[i] = NULL;
            bufsz[i] = 0;
            offt[i] = 0;
            //writers[i] = -1;
        }

        nwriters = 0;

        if(ritem->comm == gl_comps[gl_my_comp_id].comm)
            DTF_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in my comp", ritem->rank);
        else
            DTF_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in other comp", ritem->rank);


        rblock = ritem->dblocks;
        while(rblock != NULL){
            ntimes_while++;
            var_id = rblock->var_id;

            if( (witem == NULL) || (witem->var_id != var_id)){
                //find the write record for this var_id
                witem = fbuf->mst_info->iodb->witems;
                while(witem != NULL){
                    if(witem->var_id == var_id)
                        break;
                    witem = witem->next;
                }

                if(witem == NULL){
                    rblock = rblock->next;
                    /*No match right now*/
                    gl_stats.do_match_idle_time += MPI_Wtime() - t_st;
                    continue;
                }
            }

            if( (var == NULL) || (var->id != var_id))
                var = find_var(fbuf, var_id);
            assert(var != NULL);
            ndims = var->ndims;
            int nelems_to_match = 1;
            int nelems_matched;
            if(ndims > 0){
                for(i = 0; i < ndims; i++)
                    nelems_to_match *= rblock->count[i];
                matched_count = dtf_malloc(ndims*sizeof(MPI_Offset));
                assert(matched_count != NULL);
            } else
                matched_count = NULL;

            while(nelems_to_match){
                nelems_matched = 0;
                int match;
                wblock = witem->dblocks;
                if(ndims > 0)
                    while(wblock != NULL){
                        match = 0;
                        for(i = 0; i < ndims; i++)
                            if( (wblock->start[i] <= rblock->start[i]) && (rblock->start[i] < wblock->start[i] + wblock->count[i]))
                                match++;
                        if(match == ndims)
                            break;
                        wblock = wblock->next;
                    }

                if(wblock == NULL){
                    //didn't find
                    break;
                }

                n_matched_blocks++;
                //can match subblock
                matched_rank = wblock->rank;
                for(i = 0; i < ndims; i++){
                    if(rblock->start[i] + rblock->count[i] > wblock->start[i]+wblock->count[i])
                        matched_count[i] = wblock->start[i]+wblock->count[i] - rblock->start[i]; //match part
                    else
                        matched_count[i] = rblock->count[i]; //match whole
                }

                 {  /*Store*/
                    /*Find send buffer for this rank*/
                    for(i = 0; i < nwriters; i++)
                        if(writers[i] == matched_rank)
                            break;

                    if(i == nwriters){
                        /*add new send buffer*/
                        if(nwriters == allocd_nwriters ){
                            //extend buffers
                            void *tmp;
                            unsigned char **tmp1;
                            tmp = realloc((void*)writers, (nwriters+mlc_ranks)*sizeof(int));
                            assert(tmp != NULL);
                            writers = (int*)tmp;
                            gl_stats.malloc_size += mlc_ranks*sizeof(int);

                            tmp = realloc((void*)offt, (nwriters+mlc_ranks)*sizeof(size_t));
                            assert(tmp != NULL);
                            offt = (size_t*)tmp;
                            gl_stats.malloc_size += mlc_ranks*sizeof(size_t);

                            tmp = realloc((void*)bufsz, (nwriters+mlc_ranks)*sizeof(int));
                            assert(tmp != NULL);
                            bufsz = (int*)tmp;
                            gl_stats.malloc_size += mlc_ranks*sizeof(int);

                            tmp1 = realloc(sbuf, (nwriters+mlc_ranks)*sizeof(unsigned char*));
                            assert(tmp1 != NULL);
                            sbuf = tmp1;
                            gl_stats.malloc_size += mlc_ranks*sizeof(unsigned char*);
                            allocd_nwriters += mlc_ranks;

                            int j;
                            for(j = nwriters; j < allocd_nwriters; j++){
                                sbuf[j] = NULL;
                                offt[j] = 0;
                                bufsz[j] = 0;
                            }
                        }
                        offt[nwriters] = 0;
                        bufsz[nwriters] = 0;
                        sbuf[nwriters] = NULL;
                        rank_idx = nwriters;
                        writers[rank_idx] = matched_rank;
                        nwriters++;
                    } else
                        rank_idx = i;

                    /*If first time to write then, first, save the file name
                      and intracomp flag and reader rank*/
                    if(offt[rank_idx] == 0) {
                    if(bufsz[rank_idx] < MAX_FILE_NAME + sizeof(MPI_Offset)*2){
                            //extend
                            unsigned char *tmp;
                            tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + MAX_FILE_NAME + mlc_buf);
                            assert(tmp != NULL);
                            sbuf[rank_idx] = tmp;
                            bufsz[rank_idx] += MAX_FILE_NAME + mlc_buf;
                        }

                        *(MPI_Offset*)(sbuf[rank_idx]) = (MPI_Offset)ritem->rank; //to whom writer should send the data
                        offt[rank_idx] += sizeof(MPI_Offset);
                        if(intracomp_io_flag)
                            DTF_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to another w");
                        else
                           DTF_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to reader");
                        *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)intracomp_io_flag; //intra- or inter- component?
                        offt[rank_idx] += sizeof(MPI_Offset);
                        memcpy(sbuf[rank_idx] + offt[rank_idx], fbuf->file_path, MAX_FILE_NAME);//data from what file
                        offt[rank_idx] += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
                     }
                     /*Save info about the data block*/
                    //DTF_DBG(VERBOSE_ALL_LEVEL, "rank idx %d, nranks %d", rank_idx, nwriters);
                    if(offt[rank_idx] + sizeof(MPI_Offset)*ndims*2 +sizeof(MPI_Offset)> bufsz[rank_idx]){
                        unsigned char *tmp;
                        tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
                        assert(tmp != NULL);
                        sbuf[rank_idx] = tmp;
                        bufsz[rank_idx] += mlc_buf;
                    }

                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)rblock->var_id;
                    offt[rank_idx] += sizeof(MPI_Offset);
                    memcpy(sbuf[rank_idx]+offt[rank_idx], rblock->start, ndims*sizeof(MPI_Offset));
                    offt[rank_idx] += ndims*sizeof(MPI_Offset);
                    memcpy(sbuf[rank_idx]+offt[rank_idx], matched_count, ndims*sizeof(MPI_Offset));
                    offt[rank_idx] += ndims*sizeof(MPI_Offset);
                } /*store*/

                //adjust start[] and count[] in read request
                nelems_matched = 1;
                for(i = 0; i < ndims; i++){
                    rblock->start[i] += matched_count[i];
                    rblock->count[i] -= matched_count[i];
                    nelems_matched *= matched_count[i];
                }
                nelems_to_match -= nelems_matched;
            }

            dtf_free(matched_count, ndims*sizeof(MPI_Offset));

//            if(nelems_matched == 0){
//
//                rblock = rblock->next;
//                continue;
//            }

            if(nelems_to_match == 0){
                /*matched all, delete this dblock*/
                read_dblock_t *tmp = rblock;

                dtf_free(rblock->start, ndims*sizeof(MPI_Offset));
                dtf_free(rblock->count, ndims*sizeof(MPI_Offset));
                if(rblock == ritem->last_block)
                    ritem->last_block = ritem->last_block->prev;
                if(rblock == ritem->dblocks)
                    ritem->dblocks = ritem->dblocks->next;
                if(rblock->next != NULL)
                    rblock->next->prev = rblock->prev;
                if(rblock->prev != NULL)
                    rblock->prev->next = rblock->next;

                rblock = rblock->next;

                dtf_free(tmp, sizeof(read_dblock_t));
                ritem->nblocks--;
                DTF_DBG(VERBOSE_ALL_LEVEL, "Matched all in block (left %lld blocks)", ritem->nblocks);
                continue;
            }

            if(nelems_matched == 0)
                gl_stats.do_match_idle_time += MPI_Wtime() - t_st;

            rblock = rblock->next;
        }

        /*Ask writers to send the data*/
        if(nwriters > 0){
            int err;
            int my_idx = -1;

            //t_start_send = MPI_Wtime();
            for(i = 0; i < nwriters; i++){
                assert(offt[i] > 0);
                if(writers[i] == gl_my_rank){
                    my_idx = i;
                    continue;
                } else {
                    dtf_msg_t *msg = new_dtf_msg(sbuf[i], bufsz[i], IO_DATA_REQ_TAG);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Send data req to wrt %d", writers[i]);
                    err = MPI_Isend((void*)sbuf[i], offt[i], MPI_BYTE, writers[i], IO_DATA_REQ_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
                    CHECK_MPI(err);
                    ENQUEUE_ITEM(msg, gl_msg_q);
                }
            }
            //err = MPI_Waitall(nwriters, sreqs, MPI_STATUSES_IGNORE);
            //CHECK_MPI(err);
            //gl_stats.accum_comm_time += MPI_Wtime() - t_start_send;

            if(my_idx != -1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Parse data req to myself");
                send_data_wrt2rdr(sbuf[my_idx], offt[my_idx]);
                dtf_free(sbuf[my_idx], bufsz[my_idx]);
            }
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for rank %d", ritem->rank);
        }

        /*If we matched all chunks for this rank, then delete this ritem*/
        if(ritem->nblocks == 0){
            read_db_item_t *tmp = ritem;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched all. Delete ritem of rank %d (left ritems %d). ", ritem->rank,  (int)(fbuf->mst_info->iodb->nritems - 1));
            if(ritem == fbuf->mst_info->iodb->ritems)
                fbuf->mst_info->iodb->ritems = ritem->next;//fbuf->mst_info->iodb->ritems->next;
            if(ritem->prev != NULL)
                ritem->prev->next = ritem->next;
            if(ritem->next != NULL)
                ritem->next->prev = ritem->prev;
            ritem = ritem->next;
            dtf_free(tmp, sizeof(read_db_item_t));
            fbuf->mst_info->iodb->nritems--;
        } else {
//            DTF_DBG(VERBOSE_DBG_LEVEL, "Not all chunks for rreq from rank %d have been matched (%d left)", ritem->rank, (int)ritem->nchunks);
           DTF_DBG(VERBOSE_DBG_LEVEL, "Haven't matched all blocks for rank %d, %d left", ritem->rank, (int)ritem->nblocks);
           // print_read_dbitem(ritem);

            ritem = ritem->next;
        }
    }
    /*dealloc stuff*/
    dtf_free(writers, allocd_nwriters*sizeof(int));
    dtf_free(offt, allocd_nwriters*sizeof(size_t));
    dtf_free(bufsz, allocd_nwriters*sizeof(int));
    dtf_free(sbuf, allocd_nwriters*sizeof(unsigned char*));
    gl_stats.ndb_match++;

    DTF_DBG(VERBOSE_DBG_LEVEL, "after matching: %d ritems", (int)fbuf->mst_info->iodb->nritems);
}

void clean_iodb(ioreq_db_t *iodb)
{
    write_db_item_t *witem;
    read_db_item_t *ritem;
    unsigned nitems = 0, ndblocks = 0;

    witem = iodb->witems;
    while(witem != NULL){
        nitems++;

        if(witem->dblocks != NULL){
            write_dblock_t *block = witem->dblocks;
            while(block != NULL){
                ndblocks++;
                dtf_free(block->start, witem->ndims*sizeof(MPI_Offset));
                dtf_free(block->count, witem->ndims*sizeof(MPI_Offset));
                witem->dblocks = witem->dblocks->next;
                dtf_free(block, sizeof(write_dblock_t));
                witem->nblocks--;
                block = witem->dblocks;
            }
            assert(witem->nblocks == 0);
        }

        iodb->witems = iodb->witems->next;
        dtf_free(witem, sizeof(write_db_item_t));
        witem = iodb->witems;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "STATS: %u witems, %u dblocks", nitems, ndblocks);
    nitems = (unsigned)iodb->nritems;
    ndblocks = 0;
    ritem = iodb->ritems;
    while(ritem != NULL){
        if(ritem->dblocks != NULL){
            read_dblock_t *block = ritem->dblocks;
            while(block != NULL){
                ndblocks++;
                dtf_free(block->start, block->ndims*sizeof(MPI_Offset));
                dtf_free(block->count, block->ndims*sizeof(MPI_Offset));
                ritem->dblocks = ritem->dblocks->next;
                dtf_free(block, sizeof(read_dblock_t));
                block = ritem->dblocks;
                ritem->nblocks--;
            }
            assert(ritem->nblocks == 0);
        }

        iodb->ritems = iodb->ritems->next;
        dtf_free(ritem, sizeof(read_db_item_t));
        ritem = iodb->ritems;
        iodb->nritems--;
    }
//    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STATS: %u ritems, %u dblocks", nitems, ndblocks);
    iodb->witems = NULL;
    iodb->ritems = NULL;
    assert(iodb->nritems == 0);
}

static void parse_ioreqs(void *buf, int bufsz, int rank, MPI_Comm comm)
{
    int var_id, rw_flag;
    dtf_var_t *var = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    char filename[MAX_FILE_NAME];
    unsigned char *chbuf = (unsigned char*)buf;

    memcpy(filename, chbuf, MAX_FILE_NAME);
    offt += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start parsing reqs for file %s", fbuf->file_path);
    if(comm == gl_comps[fbuf->reader_id].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Reqs are from reader");
    else {
        assert(comm == gl_comps[fbuf->writer_id].comm);
        DTF_DBG(VERBOSE_DBG_LEVEL, "Req are from writer");
    }

    while(offt != (size_t)bufsz){
        rw_flag = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);

        if((var == NULL) || (var->id != var_id)){
            var = find_var(fbuf, var_id);
            assert(var != NULL);
        }

         if(rw_flag == DTF_READ){
            read_db_item_t *dbitem;
            /*Find corresponding record in the database*/
            dbitem = fbuf->mst_info->iodb->ritems;
            while(dbitem != NULL){
                if( (dbitem->rank == rank) && (dbitem->comm == comm) )
                    break;
                dbitem = dbitem->next;
            }

            if(dbitem == NULL){
                dbitem = (read_db_item_t*)dtf_malloc(sizeof(read_db_item_t));
                assert(dbitem != NULL);
                dbitem->rank = rank;
                dbitem->comm = comm;
                dbitem->next = NULL;
                dbitem->prev = NULL;
                dbitem->dblocks = NULL;
                dbitem->nblocks = 0;
                dbitem->last_block = NULL;
                //enqueue
                if(fbuf->mst_info->iodb->ritems == NULL)
                    fbuf->mst_info->iodb->ritems = dbitem;
                else{
                    dbitem->next = fbuf->mst_info->iodb->ritems;
                    fbuf->mst_info->iodb->ritems->prev = dbitem;
                    fbuf->mst_info->iodb->ritems = dbitem;
                }
                fbuf->mst_info->iodb->nritems++;
            }

            read_dblock_t *dblock = dtf_malloc(sizeof(read_dblock_t));
            assert(dblock != NULL);
            dblock->var_id = var_id;
            dblock->ndims = var->ndims;
            dblock->next = NULL;
            dblock->prev = NULL;
            dblock->start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
            assert(dblock->start != NULL);
            dblock->count = dtf_malloc(var->ndims*sizeof(MPI_Offset));
            assert(dblock->count != NULL);
            memcpy(dblock->start, chbuf+offt, var->ndims*sizeof(MPI_Offset));
            offt += var->ndims*sizeof(MPI_Offset);
            memcpy(dblock->count, chbuf+offt, var->ndims*sizeof(MPI_Offset));
            offt += var->ndims*sizeof(MPI_Offset);

            /*add to list*/
            if(dbitem->dblocks == NULL){
                dbitem->dblocks = dblock;
                dbitem->last_block = dblock;
            } else {
                dbitem->last_block->next = dblock;
                dblock->prev = dbitem->last_block;
                dbitem->last_block = dblock;
            }
            dbitem->nblocks++;
            DTF_DBG(VERBOSE_DBG_LEVEL, "ritems %d, cur item %lld blocks",(int)fbuf->mst_info->iodb->nritems, dbitem->nblocks);
        } else { /*DTF_WRITE*/
            write_db_item_t *dbitem;
            /*Allow write requests only from the writer*/
            assert(comm == gl_comps[fbuf->writer_id].comm);

            /*Find corresponding record in the database*/
            assert(fbuf->mst_info->iodb!=NULL);
            dbitem = fbuf->mst_info->iodb->witems;
            while(dbitem != NULL){
                if(dbitem->var_id == var_id)
                    break;
                dbitem = dbitem->next;
            }

            if(dbitem == NULL){
                dbitem = (write_db_item_t*)dtf_malloc(sizeof(write_db_item_t));
                assert(dbitem != NULL);
                dbitem->var_id = var_id;
                dbitem->next = NULL;
                dbitem->dblocks = NULL;
                dbitem->last_block = NULL;
                dbitem->nblocks = 0;
                dbitem->ndims = var->ndims;
                //enqueue
                if(fbuf->mst_info->iodb->witems == NULL)
                    fbuf->mst_info->iodb->witems = dbitem;
                else{
                    dbitem->next = fbuf->mst_info->iodb->witems;
                    fbuf->mst_info->iodb->witems = dbitem;
                }
            }
            write_dblock_t *dblock = dtf_malloc(sizeof(write_dblock_t));
            assert(dblock != NULL);
            dblock->rank = rank;
            dblock->next = NULL;
            if(var->ndims == 0){
                dblock->start = NULL;
                dblock->count = NULL;
            } else {
                dblock->start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                assert(dblock->start != NULL);
                dblock->count = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                assert(dblock->count != NULL);
                memcpy(dblock->start, chbuf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);
                memcpy(dblock->count, chbuf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);
            }

            /*add to list*/
            if(dbitem->dblocks == NULL){
                dbitem->dblocks = dblock;
                dbitem->last_block = dblock;
            } else {
                if(gl_conf.detect_overlap_flag){
                    /*detect if there is write overlap
                    among different ranks*/
                    write_dblock_t *tmp = dbitem->dblocks;
                    int overlap;
                    int i;
                    while(tmp != NULL){
                        overlap = 0;
                        for(i = 0; i < var->ndims; i++)
                            if( (dblock->start[i] >= tmp->start[i]) && (dblock->start[i] < tmp->start[i] + tmp->count[i]))
                                overlap++;
                            else
                                break;
                        if( (overlap == var->ndims) && (dblock->rank != tmp->rank) && (var->ndims > 0)){
                            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: data for var %d written by ranks %d \
                            and %d overlaps. State of file %s is undefined. Overlapping blocks:", var->id, dblock->rank, tmp->rank, fbuf->file_path );
                            for(i = 0; i < var->ndims; i++)
                                DTF_DBG(VERBOSE_DBG_LEVEL, "(%lld, %lld) and (%lld, %lld)",
                                dblock->start[i], dblock->count[i], tmp->start[i], tmp->count[i]);
                        }
                        tmp = tmp->next;
                    }
                }
                dbitem->last_block->next = dblock;
                dbitem->last_block = dblock;
            }
            dbitem->nblocks++;
        }
    }
    assert(offt == (size_t)bufsz);
    fbuf->mst_info->iodb->updated_flag = 1;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished parsing reqs. (mem %lu)", gl_stats.malloc_size);
}

io_req_t *new_ioreq(int id,
                    int var_id,
                    int ndims,
                    MPI_Datatype dtype,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    void *buf,
                    int rw_flag,
                    int buffered)
{

    io_req_t *ioreq = (io_req_t*)dtf_malloc(sizeof(io_req_t));
    assert(ioreq != NULL);
    int el_sz;
    MPI_Type_size(dtype, &el_sz);

    if(ndims > 0){
        int i;
        ioreq->user_buf_sz = count[0];
        for(i=1;i<ndims;i++)
            ioreq->user_buf_sz *= count[i];
//        char *s = getenv("BT_BENCH");
//        int bt_bench = 0;
//        if(s!=NULL)
//            bt_bench = 1;
//        if(bt_bench){
//            int def_el_sz;
//            MPI_Type_size(MPI_DOUBLE, &def_el_sz);
//            ioreq->user_buf_sz *= def_el_sz;
//        } else
            ioreq->user_buf_sz *= el_sz;

        ioreq->start = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->start != NULL);
        memcpy((void*)ioreq->start, (void*)start, sizeof(MPI_Offset)*ndims);
        ioreq->count = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->count != NULL);
        memcpy((void*)ioreq->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else{
        ioreq->start = NULL;
        ioreq->count = NULL;
        ioreq->user_buf_sz = el_sz;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, user bufsz %d", id, var_id, (int)ioreq->user_buf_sz);
    ioreq->is_buffered = buffered;
    if(buffered && (rw_flag == DTF_WRITE)){
        double t_start = MPI_Wtime();
        ioreq->user_buf = dtf_malloc((size_t)ioreq->user_buf_sz);
        assert(ioreq->user_buf != NULL);
        memcpy(ioreq->user_buf, buf, (size_t)ioreq->user_buf_sz);
        gl_stats.accum_dbuff_time += MPI_Wtime() - t_start;
        gl_stats.accum_dbuff_sz += (size_t)ioreq->user_buf_sz;
    } else
        ioreq->user_buf = buf;
    ioreq->next = NULL;
    ioreq->sent_flag = 0;
    ioreq->id = id;
    ioreq->var_id = var_id;
    ioreq->get_sz = 0;
    ioreq->prev = NULL;
    ioreq->dtype = dtype;
    ioreq->rw_flag = rw_flag;

    if( (rw_flag == DTF_WRITE) && gl_conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
        ioreq->checksum = compute_checksum(buf, ndims, count, dtype);
        DTF_DBG(VERBOSE_DBG_LEVEL, "checksum %.4f", ioreq->checksum);
    } else
        ioreq->checksum = 0;
    gl_stats.nioreqs++;
    return ioreq;
}

void delete_ioreq(file_buffer_t *fbuf, io_req_t **ioreq)
{
    dtf_var_t *var = find_var(fbuf, (*ioreq)->var_id);
    assert(var != NULL);

    if( (*ioreq)->is_buffered)
        dtf_free((*ioreq)->user_buf, (size_t) (*ioreq)->user_buf_sz);

    if(*ioreq == fbuf->ioreqs)
        fbuf->ioreqs = fbuf->ioreqs->next;
    if((*ioreq)->next != NULL)
        (*ioreq)->next->prev = (*ioreq)->prev;
    if( (*ioreq)->prev != NULL)
        (*ioreq)->prev->next = (*ioreq)->next;

    if( (*ioreq)->rw_flag == DTF_READ)
        fbuf->rreq_cnt--;
    else
        fbuf->wreq_cnt--;

    if( (*ioreq)->count != NULL )
        dtf_free((*ioreq)->count, var->ndims*sizeof(MPI_Offset));
    if((*ioreq)->start != NULL)
        dtf_free((*ioreq)->start, var->ndims*sizeof(MPI_Offset));
    dtf_free((*ioreq), sizeof(io_req_t));
}

/*This version divides the variable data among masters based on its
  first dimention */
void send_ioreqs_ver2(file_buffer_t *fbuf, int intracomp_match)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int mst = 0;    //only one master for now
    int nrreqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    size_t mlc_chunk = 512*1024;
    int data_to_send = 0;
    int nmasters = fbuf->mst_info->nmasters;
    int idx, err, i;
    MPI_Offset block_range, shift;
    MPI_Offset *strt, *cnt;

    if(fbuf->ioreqs == NULL)
        return;
    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    assert(sbuf != NULL);
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(offt != NULL);
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(bufsz != NULL);
    /*Distribute ioreqs between the masters based on the first dimension of the var*/

    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        bufsz[mst] = 0;
        offt[mst] = 0;
        sbuf[mst] = NULL;
    }

    nrreqs = 0;
    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            //All the following reqs should have been sent already
            break;
        }
        data_to_send = 1;

        if( (fbuf->writer_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ) && intracomp_match)
            nrreqs++;
        else if((fbuf->reader_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ))
            nrreqs++;
        ioreq = ioreq->next;
    }

    if(nrreqs == 0){
        if((gl_my_comp_id == fbuf->writer_id) && intracomp_match){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
            dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
            int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_msg_q);
        }
         /*writer will set this flag when the master root says so*/
        if(gl_my_comp_id == fbuf->reader_id)
            fbuf->done_matching_flag = 1;
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            break;
        }
        if( (var == NULL) || (var->id != ioreq->var_id)){
            var = find_var(fbuf, ioreq->var_id);
            assert(var != NULL);
        }

        DTF_DBG(VERBOSE_DBG_LEVEL, "Will send ioreq:");
        for(i = 0; i < var->ndims; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "%lld  ->  %lld", ioreq->start[i], ioreq->count[i]);

        if(var->ndims == 0){
            mst = ioreq->var_id % fbuf->mst_info->nmasters;
            if(bufsz[mst] == 0){
                    sbuf[mst] = dtf_malloc(MAX_FILE_NAME + mlc_chunk);
                    assert(sbuf[mst] != NULL);
                    bufsz[mst] = MAX_FILE_NAME + mlc_chunk;

                    memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
                    offt[mst] = MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
            }
            /*Only save the rw flag and the var id*/
            *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
            offt[mst] += sizeof(MPI_Offset);
            *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->var_id;
            offt[mst] += sizeof(MPI_Offset);
        } else {
            if(gl_conf.block_sz_range == AUTO_BLOCK_SZ_RANGE){ //divide equally among masters by the 1st dimension
                //If the first dimension is unlimited then will divide by default blocks
                if(var->shape[0] == DTF_UNLIMITED)
                    block_range = DEFAULT_BLOCK_SZ_RANGE;
                else{
                    block_range = (MPI_Offset)(var->shape[0]/fbuf->mst_info->nmasters);
                    if(block_range == 0)
                        block_range = var->shape[0];
                }
            } else
                block_range = gl_conf.block_sz_range;

            shift = 0;
            while(ioreq->start[0] + shift < ioreq->start[0] + ioreq->count[0]){
                mst = (MPI_Offset)((ioreq->start[0] + shift)/block_range)%fbuf->mst_info->nmasters;
                DTF_DBG(VERBOSE_DBG_LEVEL, "mst %d", mst);
                if(bufsz[mst] == 0){
                    sbuf[mst] = dtf_malloc(MAX_FILE_NAME + mlc_chunk);
                    assert(sbuf[mst] != NULL);
                    bufsz[mst] = MAX_FILE_NAME + mlc_chunk;

                    memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
                    offt[mst] = MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
                }

                if(var->ndims*sizeof(MPI_Offset)*3+offt[mst] > bufsz[mst]){
                    size_t ext_sz = mlc_chunk;
                    while(bufsz[mst]+ext_sz < var->ndims*sizeof(MPI_Offset)*3+offt[mst])
                        ext_sz += mlc_chunk;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu, ext sz %lu", bufsz[mst], ext_sz );
                    void *tmp = realloc(sbuf[mst], bufsz[mst] + ext_sz);
                    assert(tmp != NULL);
                    gl_stats.malloc_size += ext_sz;
                    bufsz[mst] += ext_sz;
                }

                /*Store var_id, rw_flag, start[] and count[]*/
                *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
                offt[mst] += sizeof(MPI_Offset);
                *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->var_id;
                offt[mst] += sizeof(MPI_Offset);
                memcpy(sbuf[mst]+offt[mst], ioreq->start, var->ndims*sizeof(MPI_Offset));
                //change the coord of the first dim
                *(MPI_Offset*)(sbuf[mst]+offt[mst]) = ioreq->start[0] + shift;
                strt = (MPI_Offset*)(sbuf[mst]+offt[mst]);
                offt[mst] += var->ndims*sizeof(MPI_Offset);
                memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));
                if(ioreq->count[0] - shift > block_range )
                    *(MPI_Offset*)(sbuf[mst]+offt[mst]) = block_range;
                else
                    *(MPI_Offset*)(sbuf[mst]+offt[mst]) = ioreq->count[0] - shift;
                cnt = (MPI_Offset*)(sbuf[mst]+offt[mst]);
                offt[mst] += var->ndims*sizeof(MPI_Offset);

                DTF_DBG(VERBOSE_DBG_LEVEL, "Will send info to mst %d about block:", mst);
                for(i = 0; i < var->ndims; i++)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld  ->  %lld", strt[i], cnt[i]);
                shift += block_range;
            }
            ioreq->sent_flag = 1;
        }
        ioreq = ioreq->next;
    }
    idx = -1;

    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0){
            continue;
        }
        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[mst] == gl_my_rank)){
            idx = mst;
        } else {
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], IO_REQS_TAG);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu)", fbuf->mst_info->masters[mst], offt[mst], bufsz[mst]);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, fbuf->mst_info->masters[mst], IO_REQS_TAG,
                            gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_msg_q);
        }
    }

    //gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;

    if(idx != -1){
        parse_ioreqs(sbuf[idx], (int)offt[idx], gl_my_rank, gl_comps[gl_my_comp_id].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    }

fn_exit:
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz, nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}




void send_ioreqs(file_buffer_t *fbuf, int intracomp_match)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int mst = 0;    //only one master for now
    int nrreqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    int data_to_send = 0;
    int nmasters = fbuf->mst_info->nmasters;
    int idx, err;
    MPI_Request *reqs = dtf_malloc(nmasters * sizeof(MPI_Request));
    assert(reqs != NULL);

    if(fbuf->ioreqs == NULL)
        return;
    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    assert(sbuf != NULL);
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(offt != NULL);
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(bufsz != NULL);

    /*Distribute ioreqs between the masters based on var id*/

    //alloc mem
    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        bufsz[mst] = 0;
        sbuf[mst] = NULL;
    }

    nrreqs = 0;
    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            ioreq = ioreq->next;
            //All the following reqs should have been sent already
            break;
        }
        data_to_send = 1;

        if( (fbuf->writer_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ) && intracomp_match)
            nrreqs++;
        else if((fbuf->reader_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ))
            nrreqs++;

        if( (var == NULL) || (var->id != ioreq->var_id)){
            var = find_var(fbuf, ioreq->var_id);
            assert(var != NULL);
        }
        mst = ioreq->var_id % fbuf->mst_info->nmasters;
        bufsz[mst] += sizeof(MPI_Offset)*2 + var->ndims*2*sizeof(MPI_Offset);
        ioreq = ioreq->next;
    }
    if(nrreqs == 0){
        if((gl_my_comp_id == fbuf->writer_id) && intracomp_match){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
            dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
            int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_msg_q);
        }
         /*writer will set this flag when the master root says so*/
        if(gl_my_comp_id == fbuf->reader_id)
            fbuf->done_matching_flag = 1;
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;

        bufsz[mst] += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
        DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu for mst %d", bufsz[mst], mst);
        sbuf[mst] = dtf_malloc(bufsz[mst]);
        assert(sbuf[mst] != NULL);

        memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
        offt[mst] = MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    }

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            ioreq = ioreq->next;
            //All the following reqs should have been sent already
            break;
        }
        if( (var == NULL) || (var->id != ioreq->var_id)){
            var = find_var(fbuf, ioreq->var_id);
            assert(var != NULL);
        }
        mst = ioreq->var_id % fbuf->mst_info->nmasters;
        /*Store var_id, rw_flag, start[] and count[]*/
        *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
        offt[mst] += sizeof(MPI_Offset);
        *(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->var_id;
        offt[mst] += sizeof(MPI_Offset);
        memcpy(sbuf[mst]+offt[mst], ioreq->start, var->ndims*sizeof(MPI_Offset));
        offt[mst] += var->ndims*sizeof(MPI_Offset);
        memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));
        offt[mst] += var->ndims*sizeof(MPI_Offset);
        ioreq->sent_flag = 1;
        ioreq = ioreq->next;
    }

    idx = -1;

    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0){
            reqs[mst] = MPI_REQUEST_NULL;
            continue;
        }
        assert(offt[mst] == bufsz[mst]);

        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[mst] == gl_my_rank)){
            idx = mst;
            reqs[mst] = MPI_REQUEST_NULL;
        } else {
            //dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], IO_REQS_TAG);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu)", fbuf->mst_info->masters[mst], offt[mst]);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, fbuf->mst_info->masters[mst], IO_REQS_TAG,
                            gl_comps[fbuf->writer_id].comm, &reqs[mst]);//&(msg->req));
            CHECK_MPI(err);
            //ENQUEUE_ITEM(msg, gl_msg_q);
        }
    }
    int flag = 0;
    while(!flag){
        err = MPI_Testall(nmasters, reqs, &flag, MPI_STATUSES_IGNORE);
        CHECK_MPI(err);
        progress_io_matching();
    }
    //gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;

    if(idx != -1){
        parse_ioreqs(sbuf[idx], (int)offt[idx], gl_my_rank, gl_comps[gl_my_comp_id].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    }

fn_exit:
    dtf_free(reqs, nmasters*sizeof(MPI_Request));
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz,nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}


static void reset_fbuf_match(file_buffer_t *fbuf)
{
    if ((fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag)
        /*Reset flag for future matchings*/
        fbuf->mst_info->iodb->updated_flag = 1;

    fbuf->is_matching_flag = 0;
    fbuf->mst_info->nwranks_completed = 0;
}


void notify_complete_multiple(file_buffer_t *fbuf)
{
    /*Barrier on upper level*/
    int i;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Will notify writer masters that completed multiple for %s", fbuf->file_path);
    if(fbuf->root_reader == gl_my_rank){
        for(i = 0; i < fbuf->mst_info->nmasters; i++){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Notify mst %d", fbuf->mst_info->masters[i]);
            int err = MPI_Send(fbuf->file_path,MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->masters[i], DONE_MULTIPLE_FLAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(err);
        }
    }
}



int match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag)
{

//    double t_all = MPI_Wtime();
//    double t_part1 = MPI_Wtime();
//    double t_part2=MPI_Wtime();
//    double t_part3=MPI_Wtime();
    double t_start;
    int i;
    double t_st;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs for file %s (ncid %d), intracomp %d", fbuf->file_path, fbuf->ncid, intracomp_io_flag);

    if(fbuf->mst_info->is_master_flag){
        //there should be no unfinished matches
        assert(fbuf->mst_info->nwranks_completed == 0);
    }
    assert(!fbuf->is_matching_flag);

    fbuf->done_matching_flag = 0;
    t_start = MPI_Wtime();
    /*If a writer process doesn't have any io requests, it still has to
      wait for the master process to let it complete.
      If a reader process does not have any read requests,
      it notifies the master that it completed matching and returns.*/
    DTF_DBG(VERBOSE_DBG_LEVEL, "Total %d rreqs and %d wreqs", fbuf->rreq_cnt, fbuf->wreq_cnt);
    if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
        send_ioreqs(fbuf, intracomp_io_flag);
    else /*IODB_BUILD_RANGE*/
        send_ioreqs_ver2(fbuf, intracomp_io_flag);

    gl_stats.accum_send_ioreq_time += MPI_Wtime() - t_start;

//    t_part1 = MPI_Wtime() - t_part1;
//    t_part2 = MPI_Wtime();

    fbuf->is_matching_flag = 1;
    //counter = 0;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching phase");
    while(!fbuf->done_matching_flag){
            //counter++;
           // if(counter % 500 == 0){
                t_st = MPI_Wtime();
                progress_io_matching();
               // counter = 0;
                gl_stats.accum_progr_time += MPI_Wtime() - t_st;
                t_st = MPI_Wtime();
                if( (fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag  )
                    do_matching(fbuf, intracomp_io_flag);
                gl_stats.accum_do_matching_time += MPI_Wtime() - t_st;
            //}
    }

    MPI_Barrier(fbuf->comm);

//    t_part2 = MPI_Wtime() - t_part2;
//    t_part3 = MPI_Wtime();

    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);

    if(gl_my_comp_id == fbuf->reader_id){

        if(fbuf->root_reader == gl_my_rank){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer masters reader completed matching");
            for(i = 0; i < fbuf->mst_info->nmasters; i++){
                int err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->masters[i], READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
                CHECK_MPI(err);
            }
        }
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Stat: Time for matching %.4f", MPI_Wtime() - t_start);
    reset_fbuf_match(fbuf);
    gl_stats.accum_match_time += MPI_Wtime() - t_start;

//    t_part3 = MPI_Wtime() - t_part3;
//
//    t_all = MPI_Wtime() - t_all;
//
//    DTF_DBG(VERBOSE_DBG_LEVEL, "statsss: part1: %.4f: %.4f: p2: %.4f: %.4f:prog: %.4f: %.4f: dom: %.4f: %.4f: p3: %.4f:%.4f: all: %.4f",
//            t_part1, (t_part1/t_all)*100, t_part2, (t_part2/t_all)*100, t_progress, (t_progress/t_part2)*100, t_dom, (t_dom/t_part2)*100, t_part3, (t_part3/t_all)*100, t_all);

    return 0;
}


/*writer->reader || writer->writer*/
static void send_data_wrt2rdr(void* buf, int bufsz)
{
    int var_id, rdr_rank, err, i;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    dtf_var_t *var = NULL;
    size_t rofft = 0, sofft = 0;
    unsigned char *sbuf = NULL;
    size_t sbufsz = 0;
    int def_el_sz, intracomp_flag;
    unsigned char *rbuf = (unsigned char*)buf;
    MPI_Offset *start, *count;
    int nblocks_written = 0;
    MPI_Offset *new_count = NULL, *new_start = NULL, *tmp;
    MPI_Offset nelems, fit_nelems;
    double t_start;
    int min_bl_ndims;
    size_t min_bl_sz;
    char filename[MAX_FILE_NAME];

    rdr_rank = (int)(*(MPI_Offset*)rbuf);
    rofft += sizeof(MPI_Offset);
    intracomp_flag = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Sending data to rank %d (intracomp flag %d)", rdr_rank, intracomp_flag);
    memcpy(filename, rbuf+rofft, MAX_FILE_NAME);
    rofft += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);

    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);

    t_start = MPI_Wtime();

    assert(gl_msg_buf != NULL);
    sbuf = (unsigned char*)gl_msg_buf;
    sbufsz = (int)gl_conf.data_msg_size_limit;

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: use data buf of sz %d", (int)sbufsz);

    while(rofft != (size_t)bufsz){

        if(sofft == 0){
            memcpy(sbuf, filename, MAX_FILE_NAME);
            sofft = MAX_FILE_NAME+MAX_FILE_NAME%sizeof(MPI_Offset);
        }

        var_id = (int)(*(MPI_Offset*)(rbuf+rofft));
        rofft += sizeof(MPI_Offset);

        if((var == NULL) || (var->id != var_id)){
            var = find_var(fbuf, var_id);
            assert(var != NULL);
            MPI_Type_size(var->dtype, &def_el_sz);
            /*Make sure that we can fit in at least one
              element*/
            if( (var->ndims*sizeof(MPI_Offset)+1)*2 + def_el_sz+/*padding*/def_el_sz%sizeof(MPI_Offset) > sbufsz){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: data message upper limit is \
                        too small (current value %d). Please increase by setting up DFT_DATA_MSG_SIZE_LIMIT", gl_conf.data_msg_size_limit);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
        }

        /*Prepare/Init things*/
        start = (MPI_Offset*)(rbuf+rofft);
        rofft += var->ndims*sizeof(MPI_Offset);
        count = (MPI_Offset*)(rbuf+rofft);
        rofft += var->ndims*sizeof(MPI_Offset);
        DTF_DBG(VERBOSE_DBG_LEVEL, "var id %d. Start -> count:", var->id);
        for(i = 0; i < var->ndims; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "       %lld --> %lld", start[i], count[i]);

        /*Find the ioreq that has info about this block*/
        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if( (ioreq->var_id == var_id) && (ioreq->rw_flag == DTF_WRITE)){
                int match = 0;

                if(var->ndims == 0){
                    match = 1;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Matched ioreq for scalar var %d", var->id);
                    break;
                }

                for(i = 0; i < var->ndims; i++)
                    if( (start[i] >= ioreq->start[i]) && (start[i] < ioreq->start[i]+ioreq->count[i]))
                        match++;
                    else
                        break;

                if(match == var->ndims){

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Matched ioreq with userbuf %p (strt->cnt): ", ioreq->user_buf);
                    for(i = 0; i < var->ndims; i++){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld", ioreq->start[i], ioreq->count[i]);
                        assert(start[i] + count[i] <= ioreq->start[i] + ioreq->count[i]);

                    }
                    break;
                }
            }
            ioreq = ioreq->next;
        }
        if(ioreq == NULL){
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: Matching request not found.");
        }
        assert(ioreq != NULL);
        if(gl_conf.do_checksum){
            double chsum = compute_checksum(ioreq->user_buf, var->ndims, ioreq->count, ioreq->dtype);
            if(chsum != ioreq->checksum)
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: ioreq checksum does not match with the old value: new %.4f, old %.4f", chsum, ioreq->checksum);

        }

        if(var->ndims == 0){
            min_bl_sz = def_el_sz;
            nelems = 1;
        }else {
            size_t max_mem;
            int max_nels;
            /*For subblock that will fit in the buffer.*/
            tmp = realloc(new_count, var->ndims*sizeof(MPI_Offset)); assert(tmp != NULL);new_count = tmp;
            tmp = realloc(new_start, var->ndims*sizeof(MPI_Offset)); assert(tmp != NULL);new_start = tmp;
            for(i = 0; i < var->ndims; i++)
                new_start[i] = start[i];


            /*How many elements to copy in total*/
            nelems = count[0];
            for(i = 1; i < var->ndims; i++)
                nelems *= count[i];

            /*Define the biggest full subblock that we could fit
              if the buffer was empty*/
            max_mem = sbufsz - (MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset)) -
                             var->ndims*sizeof(MPI_Offset)*2 - sizeof(MPI_Offset);
            min_bl_ndims = var->ndims;
            while(min_bl_ndims > 0){
                max_nels = 1;
                for(i = var->ndims - min_bl_ndims; i < var->ndims; i++)
                    max_nels *= count[i];
                if(max_nels*def_el_sz <= max_mem)
                    break;
                else
                    min_bl_ndims--;
            }
            min_bl_sz = max_nels*def_el_sz;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Will fit data in min blocks of sz %ld ndims %d out of %d dims", min_bl_sz, min_bl_ndims, var->ndims);
            assert(min_bl_ndims>0);
        }

        fit_nelems = 0;

        while(nelems > 0){

            size_t left_sbufsz;

            if(sofft == 0){
                memcpy(sbuf, filename, MAX_FILE_NAME);
                sofft = MAX_FILE_NAME+MAX_FILE_NAME%sizeof(MPI_Offset);
            }

            /*How much data can we fit right now?*/
            left_sbufsz = sbufsz - sofft - var->ndims*sizeof(MPI_Offset)*2 - sizeof(MPI_Offset);

            if( left_sbufsz < min_bl_sz)
                fit_nelems = 0;
            else if(var->ndims == 0)
                fit_nelems = 1;
            else {
                //adjust new_count
                for(i = 0; i < var->ndims - min_bl_ndims; i++)
                    new_count[i] = 1;
                for(i = var->ndims - min_bl_ndims; i < var->ndims; i++)
                    new_count[i] = count[i];
                /*See if we can fit several min blocks*/
                if(var->ndims != min_bl_ndims){
                    int n_fit_bl;
                    int tmpdim = var->ndims - min_bl_ndims - 1;
                    for(n_fit_bl = 2; n_fit_bl < count[tmpdim]+1 - (new_start[tmpdim] - start[tmpdim]); n_fit_bl++){
                        if( n_fit_bl*min_bl_sz <= left_sbufsz )
                            new_count[tmpdim] = n_fit_bl;
                        else
                            break;
                    }
                }
                fit_nelems = 1;
                for(i = 0; i < var->ndims; i++)
                    fit_nelems *= new_count[i];
            }

            if(fit_nelems == 0){
                assert(sofft != MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset));
                DTF_DBG(VERBOSE_DBG_LEVEL, "Send msg to %d, size %d at %.3f", rdr_rank,(int)sofft, MPI_Wtime() - gl_stats.walltime);
                /*Send this message*/
                MPI_Request req;
                int dsz = (int)sofft;
                double t_start_comm = MPI_Wtime();
                MPI_Comm comm = intracomp_flag ? gl_comps[gl_my_comp_id].comm : gl_comps[fbuf->reader_id].comm;

                err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, comm, &req);
                CHECK_MPI(err);
                err = MPI_Wait(&req, MPI_STATUS_IGNORE);
                CHECK_MPI(err);

                gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                gl_stats.ndata_msg_sent++;
                gl_stats.data_msg_sz += sofft;

                sofft = 0;
            } else {
                int cnt_mismatch = 0;
                int type_mismatch = 0;
                //int bt_bench = 0;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Will copy subblock (strt->cnt):");
                for(i=0; i < var->ndims; i++)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld \t orig: \t %lld\t --> %lld)", new_start[i], new_count[i], start[i], count[i]);

                /*Copy the block to send buffer*/

//                    char *s = getenv("BT_BENCH");
//                    if(s != NULL)  //HACK: don't check  for type mismatch for BT benchmark
//                        bt_bench = 1;
//                    if(!bt_bench && (var->dtype != ioreq->dtype)){
                if(var->dtype != ioreq->dtype){
                    /*Will only support conversion from MPI_DOUBLE to MPI_FLOAT*/
                    if(var->dtype == MPI_FLOAT){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Convert from float to double");
                        assert(ioreq->dtype == MPI_DOUBLE);
                    } else if(var->dtype == MPI_DOUBLE){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Convert from double to float");
                        assert(ioreq->dtype == MPI_FLOAT);
                    } else {
                        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF does not support this kind of type conversion");
                        assert(0);
                    }
                    type_mismatch = 1;
                }


               /*copy the subblock to the sbuf*/
                /*save var_id, start[], count[]*/
                *(MPI_Offset*)(sbuf+sofft) = (MPI_Offset)var_id;
                sofft += sizeof(MPI_Offset);
                memcpy(sbuf+sofft, new_start, var->ndims*sizeof(MPI_Offset));
                sofft += var->ndims*sizeof(MPI_Offset);
                memcpy(sbuf+sofft, new_count, var->ndims*sizeof(MPI_Offset));
                sofft += var->ndims*sizeof(MPI_Offset);

                nblocks_written++;

                DTF_DBG(VERBOSE_DBG_LEVEL, "total data to getput %lld (of nelems %lld)", fit_nelems, nelems);

                /*Copy data*/
                if(var->ndims == 0){
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Copy scalar variable");
                    get_put_data(var, ioreq->dtype, ioreq->user_buf, NULL,
                                NULL, NULL, NULL, sbuf+sofft, DTF_READ, type_mismatch);
                } else {
                    for(i = 0; i < var->ndims; i++)
                        if( (new_count[i] != 1) && (new_count[i] < ioreq->count[i]))
                            cnt_mismatch++;

                    if(cnt_mismatch <= 1){ //it's a continious block of memory
                        int req_el_sz;
                        MPI_Offset start_cpy_offt;
                        MPI_Type_size(ioreq->dtype, &req_el_sz);
                        //TODO do not include this in git (special for bt benchmark)
//                        if(bt_bench)
//                            start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, new_start) * def_el_sz;
//                        else
                            start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, new_start) * req_el_sz;

                        if(type_mismatch){
                            convertcpy(ioreq->dtype, var->dtype, (unsigned char*)ioreq->user_buf+start_cpy_offt,
                                       (void*)(sbuf+sofft), (int)fit_nelems);
                        }else{
                            memcpy(sbuf+sofft, (unsigned char*)ioreq->user_buf+start_cpy_offt, fit_nelems*def_el_sz);
                        }
                    } else {

                        get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                           ioreq->count, new_start, new_count, sbuf+sofft, DTF_READ, type_mismatch);
                    }
                }

                gl_stats.nbl++;
                sofft += fit_nelems*def_el_sz+/*padding*/(fit_nelems*def_el_sz)%sizeof(MPI_Offset);
                nelems -= fit_nelems;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Left to get put %lld", nelems);
                assert(nelems >= 0);

                if(var->ndims > 0){
                    /*shift new_start to a new position*/
                    shift_coord(var->ndims, start, count, new_start, new_count, fit_nelems);
                    if(new_start[0] == start[0]+count[0]){
                        if(nelems != 0)
                            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: nelems not zzzzero");
                        assert(nelems == 0); //we should finish now
                    }
                }
            } /*copy subblock to sbuf*/
        } /*while(nelems>0)*/

        DTF_DBG(VERBOSE_DBG_LEVEL, "Finished copying block");
        DTF_DBG(VERBOSE_DBG_LEVEL, "rofft %d (bufsz %d), Sofft %d", (int)rofft,bufsz, (int)sofft);
        //assert(sofft == sbufsz);

    } /*while(rofft != sbufsz)*/

    if(sofft > sizeof(MPI_Offset)){
        /*Send the last message*/
        DTF_DBG(VERBOSE_DBG_LEVEL, "Send (last) msg to %d size %d", rdr_rank,(int)sofft);
        MPI_Request req;
        /* First send the size of the data, then the data itself*/
        int dsz = (int)sofft;
        double t_start_comm = MPI_Wtime();
        if(intracomp_flag){
            err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[gl_my_comp_id].comm, &req);
            CHECK_MPI(err);
        } else {
            err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm, &req);
            CHECK_MPI(err);
        }
        err = MPI_Wait(&req, MPI_STATUS_IGNORE);
        CHECK_MPI(err);
        gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
        gl_stats.ndata_msg_sent++;
        gl_stats.data_msg_sz += sofft;
    }

    free(new_start);
    free(new_count);
    gl_stats.accum_extract_data_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished sending the data");
   // dtf_free(sbuf, sbufsz); //comment out since use gl_msg_buf
}

/*writer->reader*/
static void recv_data_rdr(void* buf, int bufsz)
{
    int var_id, i, nelems;
    int def_el_sz, req_el_sz;
    dtf_var_t *var = NULL;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    int nblocks_read = 0;
    int type_mismatch;
    char filename[MAX_FILE_NAME];

    MPI_Offset *start, *count;
    void *data;

    memcpy(filename, chbuf, MAX_FILE_NAME);
    offt += MAX_FILE_NAME+MAX_FILE_NAME%sizeof(MPI_Offset);
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);

    DTF_DBG(VERBOSE_DBG_LEVEL, "Received data for file %s (ncid %d)", fbuf->file_path, fbuf->ncid);
    double t_begin = MPI_Wtime();
    while(offt != bufsz){
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        if( (var == NULL) || (var->id != var_id)){
            var = find_var(fbuf, var_id);
            assert(var != NULL);
        }
        start = (MPI_Offset*)(chbuf+offt);
        offt += var->ndims*sizeof(MPI_Offset);
        count = (MPI_Offset*)(chbuf+offt);
        offt += var->ndims*sizeof(MPI_Offset);
        data = (void*)(chbuf+offt);
        nblocks_read++;

        DTF_DBG(VERBOSE_DBG_LEVEL, "var %d, bl start->count:", var->id);
        for(i = 0; i < var->ndims; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", start[i], count[i]);

        int cnt_mismatch = 0;
        /*Find the ioreq that has info about this block*/
        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if( (ioreq->var_id == var_id) && (ioreq->rw_flag == DTF_READ)){
                int match = 0;
                if(var->ndims == 0){
                    match = 1;
                    break;
                }
                for(i = 0; i < var->ndims; i++)
                    if( (start[i] >= ioreq->start[i]) && (start[i] < ioreq->start[i]+ioreq->count[i]))
                        match++;
                    else
                        break;
                if(match == var->ndims){
                    cnt_mismatch = 0;
                    for(i = 0; i < var->ndims; i++){
                        assert(start[i]+count[i]<=ioreq->start[i]+ioreq->count[i]);
                        if(count[i] != 1 && count[i] < ioreq->count[i])
                            cnt_mismatch++;
                    }
                    break;
                }

            }
            ioreq = ioreq->next;
        }
        assert(ioreq != NULL);

        /*Copy data*/
        MPI_Type_size(ioreq->dtype, &req_el_sz);
        MPI_Type_size(var->dtype, &def_el_sz);

        type_mismatch = 0;
//        int bt_bench = 0;
//         char *s = getenv("BT_BENCH");
//        if(s != NULL)  //HACK: don't check  for type mismatch for BT benchmark
//            bt_bench = 1;
        //if(!bt_bench && (ioreq->dtype != var->dtype)){
        if(ioreq->dtype != var->dtype){
            type_mismatch = 1;
            if(var->dtype == MPI_FLOAT){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Convert from float to double");
                assert(ioreq->dtype == MPI_DOUBLE);
            } else if(var->dtype == MPI_DOUBLE){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Convert from double to float");
                assert(ioreq->dtype == MPI_FLOAT);
            } else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF does not support this kind of type conversion");
                assert(0);
            }
        }


        /*Copy data*/
        nelems = 1;
        DTF_DBG(VERBOSE_DBG_LEVEL, "Getput data:");
        for(i = 0; i < var->ndims; i++){
            nelems *= count[i];
        }

       DTF_DBG(VERBOSE_DBG_LEVEL, "Will get %d elems for var %d", nelems, var->id);
        if(var->ndims == 0){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Copy scalar variable");
            get_put_data(var, ioreq->dtype, ioreq->user_buf, NULL,
                            NULL, NULL, NULL, data, DTF_WRITE, type_mismatch);
        } else if(cnt_mismatch <= 1){
            MPI_Offset start_cpy_offt;
//                if(bt_bench)
//                    start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, start) * def_el_sz;
//                else
                start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, start) * req_el_sz;
            if(type_mismatch){
                convertcpy(var->dtype, ioreq->dtype, data, (unsigned char*)ioreq->user_buf + start_cpy_offt, nelems);
            } else
                memcpy((unsigned char*)ioreq->user_buf + start_cpy_offt, data, nelems*def_el_sz);
        } else{
            /*Copy from rbuf to user buffer*/
            get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                       ioreq->count, start, count, data, DTF_WRITE, type_mismatch);
        }
        gl_stats.nbl++;
        offt += nelems*def_el_sz +(nelems*def_el_sz)%sizeof(MPI_Offset);
//            if(bt_bench)
//                ioreq->get_sz += (MPI_Offset)(nelems*def_el_sz);
//            else
            ioreq->get_sz += (MPI_Offset)(nelems*req_el_sz);


        DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, Got %d (expect %d)", ioreq->id, ioreq->var_id, (int)ioreq->get_sz, (int)ioreq->user_buf_sz);
        assert(ioreq->get_sz<=ioreq->user_buf_sz);
        if(ioreq->get_sz == ioreq->user_buf_sz){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Complete req %d (left %d), var %d, orig start->count", ioreq->id, var->id, fbuf->rreq_cnt-1);
            for(i = 0; i < var->ndims; i++)
                DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", ioreq->start[i], ioreq->count[i]);
            if(gl_conf.do_checksum){
                double chsum = compute_checksum(ioreq->user_buf, var->ndims, ioreq->count, ioreq->dtype);
                DTF_DBG(VERBOSE_DBG_LEVEL, "Chsum for req is %.4f", chsum);
            }
            //delete this ioreq
            delete_ioreq(fbuf, &ioreq);

            if(fbuf->rreq_cnt == 0){
                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE:Completed all rreqs for file %s", fbuf->file_path);

                //Only when it's intra comp reading all writer processes have
                //to notify root writer that they finished. For reader component,
                //root reader will notify all master writers later.
                if(fbuf->writer_id == gl_my_comp_id){
                    /*Notify my master writer rank that all my read io
                    requests for this file have been completed*/
                    int err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
                    CHECK_MPI(err);
                }

                /*Process of the reader component can set this flag now.
                  Writer component will set this flag only when the master
                  rank says so.*/
                if(fbuf->reader_id == gl_my_comp_id)
                    fbuf->done_matching_flag = 1;
            }
        }
    } //while(bufsz)

    gl_stats.accum_extract_data_time += MPI_Wtime()-t_begin;
    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: time to extract the data: %.3f (%d blocks)", MPI_Wtime() - t_begin, nblocks_read);
}

/*Send the file header and info about vars to the reader when the writer finishes the def mode*/
void send_file_info(file_buffer_t *fbuf, int reader_root)
{
    void *sbuf;
    MPI_Offset sbuf_sz, err;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Will send file info to reader %d", reader_root);
    pack_file_info(fbuf, &sbuf_sz, &sbuf);
    assert(sbuf_sz > 0);
    double t_start = MPI_Wtime();
    MPI_Request req;
    err = MPI_Isend(sbuf, (int)sbuf_sz, MPI_BYTE, reader_root, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm, &req);
    CHECK_MPI(err);
    err = MPI_Wait(&req, MPI_STATUS_IGNORE);
    CHECK_MPI(err);
    gl_stats.accum_comm_time += MPI_Wtime() - t_start;
    dtf_free(sbuf, sbuf_sz);
}


static void notify_masters(file_buffer_t *fbuf, int msgtag)
{
    int i, err;
    assert(fbuf->mst_info != NULL);
    assert(gl_my_rank == fbuf->root_writer);
    for(i = 0; i < fbuf->mst_info->nmasters; i++) {
        if(fbuf->mst_info->masters[i] == gl_my_rank)
            continue;
        err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->masters[i], msgtag, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(err);
    }
}

/*Function executed by master writers to notify other writers
  about something defined by mpitag (match completion or file close)*/
static void notify_workgroup(file_buffer_t *fbuf, int msgtag)
{
    int i, err;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Mst %d will notify workgroup (msgtag %d) for %s", gl_my_rank, msgtag, fbuf->file_path);
    assert(fbuf->mst_info->is_master_flag);

    for(i = 0; i < fbuf->mst_info->my_wg_sz - 1; i++){
        err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->my_wg[i], msgtag, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(err);
    }

//    /*First, translate the ranks in the communicator which
//    was used to open the file to global ranks*/
//    MPI_Comm_group(gl_comps[gl_my_comp_id].comm, &glob_group);
//    MPI_Comm_group(fbuf->comm, &file_group);
//
//    nranks = fbuf->mst_info->my_workgroup_sz - 1;
//    ranks = dtf_malloc(nranks*sizeof(int));
//    assert(ranks != NULL);
//    glob_ranks = dtf_malloc(nranks*sizeof(int));
//    assert(glob_ranks != NULL);
//
//    MPI_Comm_rank(fbuf->comm, &rank);
//
//    for(i = 0; i < nranks; i++)
//        ranks[i] = rank + i + 1;
//
//    err = MPI_Group_translate_ranks(file_group, nranks, ranks, glob_group, glob_ranks);
//    CHECK_MPI(err);
//
//    for(i = 0; i < nranks; i++){
//        err = MPI_Send(&fbuf->ncid, 1, MPI_INT, glob_ranks[i], msgtag, gl_comps[gl_my_comp_id].comm);
//        CHECK_MPI(err);
//        gl_stats.nmsg_sent++;
//    }
//
//    MPI_Group_free(&file_group);
//    MPI_Group_free(&glob_group);
//    dtf_free(ranks, nranks*sizeof(int));
//    dtf_free(glob_ranks, nranks*sizeof(int));
}

void progress_msg_queue()
{
    dtf_msg_t *msg;
    int flag, err;
    MPI_Status status;
    double t_st;

    if(gl_msg_q == NULL)
       return;

    msg = gl_msg_q;
    while(msg != NULL){
        err = MPI_Test(&(msg->req), &flag, &status);
        CHECK_MPI(err);
        if(flag){
            t_st = MPI_Wtime();
            dtf_msg_t *tmp = msg->next;
            DEQUEUE_ITEM(msg, gl_msg_q);
            if(msg->tag == FILE_READY_TAG){
                file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, (char*)msg->buf, -1);
                if(fbuf == NULL)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Error: cant find %s", (char*)msg->buf);
                assert(fbuf != NULL);
                assert(fbuf->fready_notify_flag == RDR_NOTIF_POSTED);
                fbuf->fready_notify_flag = RDR_NOTIFIED;
            }
            delete_dtf_msg(msg);
            msg = tmp;
            gl_stats.progr_work_time += MPI_Wtime() - t_st;
        } else
            msg = msg->next;
    }
}

static void print_recv_msg(int tag)
{
    switch(tag){
        case FILE_READY_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_READY_TAG");
            break;
        case IO_DATA_REQ_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_DATA_REQ_TAG");
            break;
        case DONE_MULTIPLE_FLAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag DONE_MULTIPLE_FLAG");
            break;
        case READ_DONE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag READ_DONE_TAG");
            break;
        case IO_CLOSE_FILE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_CLOSE_FILE_TAG");
            break;
        case IO_OPEN_FILE_FLAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_OPEN_FILE_FLAG");
            break;
        case FILE_INFO_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_TAG");
            break;
        case FILE_INFO_REQ_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_REQ_TAG");
            break;
        case ROOT_MST_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag ROOT_MST_TAG");
            break;
        case MATCH_DONE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag MATCH_DONE_TAG");
            break;
        case IO_REQS_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_REQS_TAG");
            break;
        case IO_DATA_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_DATA_TAG");
            break;
        default:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag unknown %d", tag);
    }
}

void progress_io_matching()
{
    MPI_Status status;
    int comp, flag, src, err;
    file_buffer_t *fbuf;

    int bufsz;
    void *rbuf;
    char filename[MAX_FILE_NAME];
    double t_st;
    double t_start_comm;
    gl_stats.nprogress_call++;
    /*first check if there are any requests for file info
      that we can process. */
    process_file_info_req_queue();
    progress_msg_queue();

    for(comp = 0; comp < gl_ncomp; comp++){
        if( gl_comps[comp].comm == MPI_COMM_NULL){
            continue;
        }

        while(1){
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                break;
            }
            src = status.MPI_SOURCE;

            print_recv_msg(status.MPI_TAG);
            t_st = MPI_Wtime();
            switch(status.MPI_TAG){
                case FILE_INFO_REQ_TAG:
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv root req from rank %d (comp %s)", src, gl_comps[comp].name);
                    rbuf = dtf_malloc(MAX_FILE_NAME+sizeof(int));
                    assert(rbuf != NULL);
                    err = MPI_Recv(rbuf, (int)(MAX_FILE_NAME+sizeof(int)), MPI_BYTE, src, FILE_INFO_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    memcpy(filename, rbuf, MAX_FILE_NAME);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    if(fbuf == NULL)
                        DTF_DBG(VERBOSE_ERROR_LEVEL, "File not found %s", filename);
                    assert(fbuf != NULL);

                    /* If the global master zero (global rank 0) couldn't process the
                       request immediately because it doesn't know yet who is the root
                       writer for the file, the req has to be enqueue and, later, the master
                       will periodically check the queue to see if it can process the
                       request.
                     */
                    if(fbuf->root_writer == -1){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Don't know who is the root for %s now. Will queue the req.", fbuf->file_path);

                        file_info_req_q_t *req = dtf_malloc(sizeof(file_info_req_q_t));
                        assert(req != NULL);
                        memcpy(req->filename, filename, MAX_FILE_NAME);
                        req->buf = rbuf;
                        req->next = NULL;
                        req->prev = NULL;
                        ENQUEUE_ITEM(req, gl_finfo_req_q);

                        //enqueue
//                        if(gl_finfo_req_q == NULL)
//                            gl_finfo_req_q = req;
//                        else {
//                            req->next = gl_finfo_req_q;
//                            gl_finfo_req_q->prev = req;
//                            gl_finfo_req_q = req;
//                        }

                        break;
                    }
                    if(gl_my_rank == fbuf->root_writer){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");

                        if(fbuf->iomode == DTF_IO_MODE_FILE){
                            memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
                            /*just notify the reader that I am the root writer*/

                            err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm);
                            CHECK_MPI(err);
                            fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;

                        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
                            memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
                            send_file_info(fbuf, fbuf->root_reader);
                        }
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Root reader for file %s is %d", fbuf->file_path, fbuf->root_reader);
                        dtf_free(rbuf, MAX_FILE_NAME+sizeof(int));
                    } else {
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Forward the request to rank %d", fbuf->root_writer);
                        dtf_msg_t *msg = new_dtf_msg(rbuf, MAX_FILE_NAME+sizeof(int), FILE_INFO_REQ_TAG);

                        err = MPI_Isend(msg->buf,(int)(MAX_FILE_NAME+sizeof(int)), MPI_BYTE, fbuf->root_writer,
                                           FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
                        CHECK_MPI(err);
                        ENQUEUE_ITEM(msg, gl_msg_q);
                    }

                    break;
                case IO_REQS_TAG:
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Received reqs from %d, comp %d (my comp %d), bufsz %d", src, comp,gl_my_comp_id, (int)bufsz);
                    rbuf = dtf_malloc(bufsz);
                    assert(rbuf != NULL);
                    t_start_comm = MPI_Wtime();
                    err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_REQS_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    parse_ioreqs(rbuf, bufsz, src, gl_comps[comp].comm);
                    dtf_free(rbuf, bufsz);
                    break;
                case IO_DATA_REQ_TAG:
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    rbuf = dtf_malloc(bufsz);
                    assert(rbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv data req from mst %d", src);
                    t_start_comm = MPI_Wtime();
                    err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    send_data_wrt2rdr(rbuf, bufsz);
                    dtf_free(rbuf, bufsz);
                    break;

                case IO_DATA_TAG:
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recved data from %d", src);
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    assert(bufsz>0);
                    assert(gl_msg_buf != NULL);
                    rbuf = gl_msg_buf; //dtf_malloc(bufsz);
                    assert(rbuf != NULL);
                    t_start_comm = MPI_Wtime();
                    err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: recv data from %d", src);
                    recv_data_rdr(rbuf, bufsz);
                   // dtf_free(rbuf, bufsz);
                    break;
                case READ_DONE_TAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, READ_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv read done for file %s", fbuf->file_path);
                    if(comp == gl_my_comp_id){
                        assert(fbuf->root_writer == gl_my_rank);
                        fbuf->mst_info->nwranks_completed++;
                        DTF_DBG(VERBOSE_DBG_LEVEL, " -> from writer %d (tot %d)", src, fbuf->mst_info->nwranks_completed);
                        if((fbuf->mst_info->nwranks_opened > 0) && (fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened)){

                            assert(gl_my_rank == fbuf->root_writer);
                            DTF_DBG(VERBOSE_DBG_LEVEL, "Notify masters & workgroup that intracomp matching is finished");
                            notify_masters(fbuf, MATCH_DONE_TAG);
                            notify_workgroup(fbuf, MATCH_DONE_TAG);
                            fbuf->done_matching_flag = 1;
                            DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE, Done matching at %.3f", MPI_Wtime()-gl_stats.walltime);
                        }
                    } else {
                        assert(comp == fbuf->reader_id);
                        //fbuf->mst_info->nrranks_completed++;
                        DTF_DBG(VERBOSE_DBG_LEVEL, "->from reader root. Match complete. Notify my group.");
                        assert(fbuf->mst_info->nwranks_completed == 0);
                        //notify_masters(fbuf, MATCH_DONE_TAG);
                        notify_workgroup(fbuf, MATCH_DONE_TAG);
                        fbuf->done_matching_flag = 1;
                    }
                    break;
                case MATCH_DONE_TAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, MATCH_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    if(fbuf->mst_info->is_master_flag){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                        /*Tell other writer ranks that they can complete matching*/
                        notify_workgroup(fbuf, MATCH_DONE_TAG);
                    }
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                    fbuf->done_matching_flag = 1;
/*
//                    if(fbuf->rdr_closed_flag){
////                        DTF_DBG(VERBOSE_DBG_LEVEL, "Cleaning up all reqs and dbs");
//                        assert(fbuf->rreq_cnt == 0);
//
//                        delete_ioreqs(fbuf);
//
//                        if(fbuf->mst_info->is_master_flag){
//                            //Check that I don't have any read reqs incompleted
//                            assert(fbuf->mst_info->iodb->nritems == 0);
//                            //Clean my iodb
//                            clean_iodb(fbuf->mst_info->iodb);
//                        }
                    //}
                    */
                    break;
                case IO_CLOSE_FILE_TAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, IO_CLOSE_FILE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv close file tag for %s from %d", fbuf->file_path, src);
                    //assert(!fbuf->rdr_closed_flag); //check that haven't received that before

                    if(gl_my_rank == fbuf->root_writer){
                        //DTF_DBG(VERBOSE_DBG_LEVEL, "Notify other masters that readers are closing the file");
                        //notify_masters(fbuf, IO_CLOSE_FILE_TAG);

                        //TODO for multi-cycle version will need to figure out with these states and notifications
                        if(fbuf->iomode == DTF_IO_MODE_FILE)
                            fbuf->fready_notify_flag = RDR_HASNT_OPENED;
                    }
                    if(fbuf->mst_info->is_master_flag){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that they can close the file %s", fbuf->file_path);
                        notify_workgroup(fbuf, IO_CLOSE_FILE_TAG);
                    }
                    //TODO figure out with letkf multiple closing file when to reset this flag?
                    fbuf->rdr_closed_flag = 1;

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Close flag set for file %s", fbuf->file_path);
//                    if(fbuf->done_matching_flag){
//                        if(fbuf->mst_info->is_master_flag){
//                            /*Check that I don't have any read reqs incompleted*/
//                            assert(fbuf->mst_info->iodb->nritems == 0);
//                            /*Clean my iodb*/
//                            clean_iodb(fbuf->mst_info->iodb);
//                        }
//                         /*Can delete write requests only if all ps-s have finished
//                         matching.*/
//                        delete_ioreqs(fbuf);
//                    }
                    break;
                case ROOT_MST_TAG:
                    src = status.MPI_SOURCE;
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, ROOT_MST_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive ROOT_MST_TAG notif for %s", filename);

                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    assert(fbuf->root_writer == -1);
                    fbuf->root_writer = src;
                    break;
                case FILE_READY_TAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);

                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    fbuf->is_ready = 1;
                    break;
                case DONE_MULTIPLE_FLAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, DONE_MULTIPLE_FLAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv done multiple tag for %s from %d", fbuf->file_path, src);

                    if(fbuf->mst_info->is_master_flag){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that multiple matching for %s completed", fbuf->file_path);
                        /*Tell other writer ranks that they can complete matching*/
                        notify_workgroup(fbuf, DONE_MULTIPLE_FLAG);
                    }
                    fbuf->done_matching_flag = 1;
                    fbuf->done_match_multiple_flag = 1;
                    break;
                case IO_OPEN_FILE_FLAG:
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, IO_OPEN_FILE_FLAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv open file tag for %s from %d", fbuf->file_path, src);
                    if(fbuf->mst_info->is_master_flag){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify workgroup");
                        notify_workgroup(fbuf, IO_OPEN_FILE_FLAG);
                        fbuf->rdr_closed_flag = 0;

                        if(fbuf->iomode == DTF_IO_MODE_FILE)
                            fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;

                    }
                    break;
                default:
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: unknown tag %d", status.MPI_TAG);
                    assert(0);
            }
            gl_stats.progr_work_time += MPI_Wtime() - t_st;
        }
    }
}


/*function called by the writer processes*/
int init_req_match_masters(MPI_Comm comm, master_info_t *mst_info)
{
    int wg, nranks, myrank, i, err;
    char* s = getenv("MAX_WORKGROUP_SIZE");
    MPI_Group glob_group, file_group;
    int my_master, my_master_glob;
    int *masters, *ranks;

    if(s == NULL)
        wg = MAX_WORKGROUP_SIZE;
    else
        wg = atoi(s);
    assert(wg > 0);

    /* First, find out my master and, if I am a master, find out the size of my
       workgroup. Create the list of masters. */
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &myrank);

    MPI_Comm_group(gl_comps[gl_my_comp_id].comm, &glob_group);
    MPI_Comm_group(comm, &file_group);

    if(nranks <= wg){
        my_master = 0;
        mst_info->my_wg_sz = nranks;
        mst_info->nmasters = 1;
    } else {
        my_master = (int)(myrank/wg) * wg;
        mst_info->my_wg_sz = wg;
        mst_info->nmasters = (int)(nranks/wg);
        if(nranks % wg > 0){
            mst_info->nmasters++;
            if(myrank >= (mst_info->nmasters-1)*wg)
                mst_info->my_wg_sz = nranks % wg;
        }
    }

    if(myrank == 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Nmasters %d", mst_info->nmasters);

    /*Translate the rank of my master*/
    err = MPI_Group_translate_ranks(file_group, 1, &my_master, glob_group, &my_master_glob);
    CHECK_MPI(err);
    mst_info->is_master_flag = (gl_my_rank == my_master_glob) ? 1 : 0;

    mst_info->masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    assert(mst_info->masters != NULL);
    masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    assert(masters != NULL);

    masters[0] = 0;
    for(i = 1; i < mst_info->nmasters; i++){
        masters[i] = masters[i-1] + wg;
    }
    /*Translate rank from subcommunicator to the global rank in mpi_comm_world*/
    err = MPI_Group_translate_ranks(file_group, mst_info->nmasters, masters, glob_group, mst_info->masters);
    CHECK_MPI(err);

    if(myrank == 0){
        for(i = 0; i < mst_info->nmasters; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "Rank %d is a master", mst_info->masters[i]);
    }

    dtf_free(masters, mst_info->nmasters * sizeof(int));

    if(mst_info->is_master_flag){
        DTF_DBG(VERBOSE_DBG_LEVEL, "My wg size %d", mst_info->my_wg_sz);
        //Translate my workgroup ranks
        if(mst_info->my_wg_sz == 1){
            //master is the only rank in the wg
            mst_info->my_wg = NULL;
        } else {
            nranks = mst_info->my_wg_sz - 1;
            mst_info->my_wg = dtf_malloc(nranks*sizeof(int));
            assert(mst_info->my_wg != NULL);
            ranks = dtf_malloc(nranks*sizeof(int));
            assert(ranks != NULL);

            for(i = 0; i < nranks; i++)
                ranks[i] = myrank + i + 1;

            err = MPI_Group_translate_ranks(file_group, nranks, ranks, glob_group, mst_info->my_wg);
            CHECK_MPI(err);
            dtf_free(ranks, nranks*sizeof(int));
        }

    } else {
        mst_info->my_wg = NULL;
    }
    MPI_Group_free(&file_group);
    MPI_Group_free(&glob_group);

    mst_info->nwranks_completed = 0;
    mst_info->iodb = NULL;
    return 0;
}

