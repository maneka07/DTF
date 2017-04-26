#include "dtf_req_match.h"
#include "dtf_common.h"
#include "dtf_util.h"
#include "dtf.h"
#include <unistd.h>

file_info_req_q_t *gl_finfo_req_q = NULL;

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
       - file ncid
       - file header size
       - header
       - number of masters
       - master list
       - number of vars
       - vars
    */

    sz =   (MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset)) +
           (fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset)) +
           fbuf->mst_info->nmasters*sizeof(MPI_Offset) +
           sizeof(MPI_Offset)*4;
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
    /*ncid*/
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->ncid;
    offt += sizeof(MPI_Offset);
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
        memcpy(chbuf+offt, fbuf->mst_info->masters, fbuf->mst_info->nmasters+sizeof(int));
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
       - file ncid
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
    fbuf->ncid = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
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
    //NOTE: hidden bug: if ndims*sizeof(MPI_Offset) > 512 there will be segmentation fault
    //very unlikely to happen in real life though
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

    int n_matched_blocks = 0;

    double t_begin, t_findblock = 0, t_store_data=0;
    double t_start_send, t_send = 0;
    double t_start;

    if(!fbuf->mst_info->iodb->updated_flag) //no new info since last time matching was done, ignore
        return;
    fbuf->mst_info->iodb->updated_flag = 0; //reset

    t_start = MPI_Wtime();

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

        nwriters = 0;

        for(i = 0; i < allocd_nwriters; i++){
            //sbuf[i] = NULL;
            //bufsz[i] = 0;
            offt[i] = 0;
            //writers[i] = -1;
        }
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
                    continue;
                }
            }

            //DTF_DBG(VERBOSE_ALL_LEVEL, "write record for this var:");
            if( (var == NULL) || (var->id != var_id))
                var = find_var(fbuf, var_id);
            assert(var != NULL);
            ndims = var->ndims;
            //TODO what if a scalar var?
            assert(ndims > 0);

            int nelems_to_match = 1;
            int nelems_matched;
            for(i = 0; i < ndims; i++)
                nelems_to_match *= rblock->count[i];
            matched_count = dtf_malloc(ndims*sizeof(MPI_Offset));
            assert(matched_count != NULL);

            while(nelems_to_match){
                nelems_matched = 0;
                int match;

                t_begin = MPI_Wtime();
                wblock = witem->dblocks;
                while(wblock != NULL){
                    match = 0;
                    for(i = 0; i < ndims; i++)
                        if( (wblock->start[i] <= rblock->start[i]) && (rblock->start[i] < wblock->start[i] + wblock->count[i]))
                            match++;
                    if(match == ndims)
                        break;
                    wblock = wblock->next;
                }
                t_findblock += MPI_Wtime() - t_begin;

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
                  t_begin = MPI_Wtime();
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
                            gl_stats.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            writers = (int*)tmp;

                            tmp = realloc((void*)offt, (nwriters+mlc_ranks)*sizeof(size_t));
                            gl_stats.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            offt = (size_t*)tmp;

                            tmp = realloc((void*)bufsz, (nwriters+mlc_ranks)*sizeof(int));
                            gl_stats.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            bufsz = (int*)tmp;

                            tmp1 = realloc(sbuf, (nwriters+mlc_ranks)*sizeof(unsigned char*));
                            gl_stats.malloc_size += mlc_ranks*sizeof(unsigned char*);
                            assert(tmp1 != NULL);
                            sbuf = tmp1;
                            allocd_nwriters += mlc_ranks;

                            int j;
                            for(j = nwriters; j < allocd_nwriters; j++)
                                sbuf[j] = NULL;
                        }
                        offt[nwriters] = 0;
                        bufsz[nwriters] = 0;
                        sbuf[nwriters] = NULL;
                        rank_idx = nwriters;
                        writers[rank_idx] = matched_rank;
                        nwriters++;
                    } else
                        rank_idx = i;

                    /*If first time to write then, first, save the file ncid
                      and intracomp flag and reader rank*/
                    if(offt[rank_idx] == 0) {
                        if(bufsz[rank_idx] < sizeof(MPI_Offset)*3){
                            //extend
                            unsigned char *tmp;
                            tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
                            assert(tmp != NULL);
                            sbuf[rank_idx] = tmp;
                            bufsz[rank_idx] += mlc_buf;
                            gl_stats.malloc_size += mlc_buf;
                        }

                        *(MPI_Offset*)(sbuf[rank_idx]) = (MPI_Offset)ritem->rank; //to whom writer should send the data
                        offt[rank_idx] += sizeof(MPI_Offset);
                        if(intracomp_io_flag)
                            DTF_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to another w");
                        else
                           DTF_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to reader");
                        *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)intracomp_io_flag; //intra- or inter- component?
                        offt[rank_idx] += sizeof(MPI_Offset);
                        *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)fbuf->ncid;  //data from what file
                        offt[rank_idx] += sizeof(MPI_Offset);
                     }
                     /*Save info about the data block*/
                    //DTF_DBG(VERBOSE_ALL_LEVEL, "rank idx %d, nranks %d", rank_idx, nwriters);
                    if(offt[rank_idx] + sizeof(MPI_Offset)*ndims*2 +sizeof(MPI_Offset)> bufsz[rank_idx]){
                        unsigned char *tmp;
                        tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
                        assert(tmp != NULL);
                        sbuf[rank_idx] = tmp;
                        bufsz[rank_idx] += mlc_buf;
                        gl_stats.malloc_size += mlc_buf;
                    }

                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)rblock->var_id;
                    offt[rank_idx] += sizeof(MPI_Offset);
                    memcpy(sbuf[rank_idx]+offt[rank_idx], rblock->start, ndims*sizeof(MPI_Offset));
                    offt[rank_idx] += ndims*sizeof(MPI_Offset);
                    memcpy(sbuf[rank_idx]+offt[rank_idx], matched_count, ndims*sizeof(MPI_Offset));
                    offt[rank_idx] += ndims*sizeof(MPI_Offset);
                } /*store*/
                t_store_data += MPI_Wtime() - t_begin;
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

            if(nelems_matched == 0){

                rblock = rblock->next;
                continue;
            }

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

            rblock = rblock->next;
        }

        if(nwriters > 0){
            t_start_send = MPI_Wtime();
            int errno;
            MPI_Request *sreqs = (MPI_Request*)dtf_malloc(nwriters*sizeof(MPI_Request));
            assert(sreqs != NULL);
            int my_idx = -1;
            for(i = 0; i < nwriters; i++){
                assert(offt[i] > 0);

                if(writers[i] == gl_my_rank){
                    my_idx = i;
                    sreqs[i] = MPI_REQUEST_NULL;
                    continue;
                } else {
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Send data req to wrt %d", writers[i]);
                    t_start_send = MPI_Wtime();
                    errno = MPI_Isend((void*)sbuf[i], offt[i], MPI_BYTE, writers[i], IO_DATA_REQ_TAG, gl_comps[gl_my_comp_id].comm, &sreqs[i]);
                    CHECK_MPI(errno);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_send;
                }
                gl_stats.nmatching_msg_sent++;
                gl_stats.accum_msg_sz += offt[i];
            }

            errno = MPI_Waitall(nwriters, sreqs, MPI_STATUSES_IGNORE);
            CHECK_MPI(errno);
            gl_stats.accum_comm_time += MPI_Wtime() - t_start_send;
            if(my_idx != -1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Parse data req to myself");
                send_data_wrt2rdr(sbuf[my_idx], offt[my_idx]);
            }
            t_send += MPI_Wtime() - t_start_send;
//            completed = 0;
//            //DTF_DBG(VERBOSE_ALL_LEVEL, "Start waiting (rreq from rank %d)", ritem->rank);
//            double t_progress = 0;
//            while(!completed){
//                errno = MPI_Testall(nwriters, sreqs, &completed, MPI_STATUSES_IGNORE);
//                CHECK_MPI(errno);
//                double tmp = MPI_Wtime();
//                progress_io_matching();
//                t_progress += MPI_Wtime() - tmp;
//            }

//            gl_stats.accum_comm_time += MPI_Wtime() - t_start_send;// - t_progress;
//            t_send += MPI_Wtime() - t_start_send;


            //DTF_DBG(VERBOSE_ALL_LEVEL, "Finish waiting (rreq from rank %d)", ritem->rank);
            dtf_free(sreqs, nwriters*sizeof(MPI_Request));
//            for(i = 0; i < nwriters; i++)
//                dtf_free(sbuf[i], offt[i]);

            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for rank %d", ritem->rank);

        } //else {
//            DTF_DBG(VERBOSE_DBG_LEVEL, "Did not manage to match anything for rank %d", ritem->rank);
//        }

        /*If we matched all chunks for this rank, then delete this ritem*/
        if(ritem->nblocks == 0){
            read_db_item_t *tmp = ritem;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched all. Delete ritem of rank %d (left ritems %d). ", ritem->rank,  (int)(fbuf->mst_info->iodb->nritems - 1));
            //DTF_DBG(VERBOSE_ALL_LEVEL, "Ritems before:");
            //print_ritems(fbuf->mst_info->iodb->ritems);
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

    for(i = 0; i < allocd_nwriters; i++){
        if(sbuf[i] != NULL)
            dtf_free(sbuf[i], (size_t)bufsz[i]);
    }
    /*dealloc stuff*/
    dtf_free(writers, allocd_nwriters*sizeof(int));
    dtf_free(offt, allocd_nwriters*sizeof(size_t));
    dtf_free(bufsz, allocd_nwriters*sizeof(int));
    dtf_free(sbuf, allocd_nwriters*sizeof(unsigned char*));
    DTF_DBG(VERBOSE_DBG_LEVEL, "Time to send data in do_matching %.3f", t_send);
    gl_stats.ndb_match++;
    gl_stats.accum_db_match_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "after matching: %d ritems", (int)fbuf->mst_info->iodb->nritems);
    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: %d times while, %d blocks matched, t_store %.3f, t_find_block %.3f",
            ntimes_while, n_matched_blocks, t_store_data, t_findblock );

    //DTF_DBG(VERBOSE_DBG_LEVEL, "Stat: Time to match db reqs %.4f, time to send %.4f", MPI_Wtime() - t_start, t_send);
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
    int ncid, var_id, rw_flag;
    dtf_var_t *var = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    ncid = (int)(*((MPI_Offset*)chbuf));
    offt += sizeof(MPI_Offset);
    double t_start = MPI_Wtime();

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
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
    DTF_DBG(VERBOSE_DBG_LEVEL, "Stat: time to parse reqs %.4f", MPI_Wtime()-t_start);
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
    if(buffered){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Buffering the data");
        ioreq->user_buf = dtf_malloc((size_t)ioreq->user_buf_sz);
        assert(ioreq->user_buf != NULL);
        memcpy(ioreq->user_buf, buf, (size_t)ioreq->user_buf_sz);
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

void send_ioreqs(file_buffer_t *fbuf, int intracomp_match)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int mst = 0;    //only one master for now
    int n_intra_reqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    int nmasters = 1;
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
    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++)
        bufsz[mst] = 0;
    int data_to_send = 0;

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            ioreq = ioreq->next;
            continue;
        }
        data_to_send = 1;
        if( (fbuf->writer_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ))
            n_intra_reqs++;

        if( (var == NULL) || (var->id != ioreq->var_id)){
            var = find_var(fbuf, ioreq->var_id);
            assert(var != NULL);
        }
        mst = ioreq->var_id % fbuf->mst_info->nmasters;
        bufsz[mst] += sizeof(MPI_Offset)*2 + var->ndims*2*sizeof(MPI_Offset);
        ioreq = ioreq->next;
    }
    if(intracomp_match && (n_intra_reqs == 0)){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Writer has no read requests. Nofity master.");
        int errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
        CHECK_MPI(errno);
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send


    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;

        bufsz[mst] += sizeof(MPI_Offset); //ncid
        DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu for mst %d", bufsz[mst], mst);
        sbuf[mst] = dtf_malloc(bufsz[mst]);
        assert(sbuf[mst] != NULL);

        /*save ncid*/
        *(MPI_Offset*)sbuf[mst] = (MPI_Offset)fbuf->ncid;
        offt[mst] = sizeof(MPI_Offset);
    }

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            ioreq = ioreq->next;
            continue;
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

    for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;
        assert(offt[mst] == bufsz[mst]);

        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[mst] == gl_my_rank)){
                parse_ioreqs(sbuf[mst], (int)offt[mst], gl_my_rank, gl_comps[gl_my_comp_id].comm);
        } else {
            MPI_Request req;
            int errno;
            assert(gl_comps[fbuf->writer_id].comm != MPI_COMM_NULL);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Send reqs to mst %d (bufsz %lu)", fbuf->mst_info->masters[mst], offt[mst]);
            double t_start = MPI_Wtime();
            errno = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, fbuf->mst_info->masters[mst], IO_REQS_TAG, gl_comps[fbuf->writer_id].comm, &req);
            CHECK_MPI(errno);
            gl_stats.nmatching_msg_sent++;
            gl_stats.accum_msg_sz += bufsz[mst];
            errno = MPI_Wait(&req, MPI_STATUS_IGNORE);
            gl_stats.accum_comm_time += MPI_Wtime() - t_start;
            CHECK_MPI(errno);
        }
    }
fn_exit:
     for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;
        dtf_free(sbuf[mst], bufsz[mst]);
    }
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz, nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}


static void reset_fbuf_match(file_buffer_t *fbuf)
{
    if ((fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag)
        /*Reset flag for future matchings*/
        fbuf->mst_info->iodb->updated_flag = 1;

    fbuf->is_matching_flag = 0;
    fbuf->mst_info->nrranks_completed = 0;
    fbuf->mst_info->nwranks_completed = 0;
}


void notify_complete_multiple(file_buffer_t *fbuf)
{
    int i;
    MPI_Barrier(fbuf->comm);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Will notify writer masters that completed multiple for %s", fbuf->file_path);
    int rank;
    //TODO init root-_reader for readers
    MPI_Comm_rank(fbuf->comm, &rank);
    if(rank == 0){
        for(i = 0; i < fbuf->mst_info->nmasters; i++){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Notify mst %d", fbuf->mst_info->masters[i]);
            int errno = MPI_Send(&fbuf->ncid, 1, MPI_INT, fbuf->mst_info->masters[i], DONE_MULTIPLE_FLAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(errno);
        }
    }
}

void match_ioreqs_all(int rw_flag)
{
    file_buffer_t *fbuf;
    io_req_t *ioreq;
    int fbuf_cnt = 0, tmp_cnt;

    /*Note: rw_flag should be DTF_WRITE (checked this in higher level)*/
    assert(rw_flag == DTF_WRITE);

    DTF_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs all");

    /*Do preparations*/
    fbuf = gl_filebuf_list;
    while(fbuf != NULL){
        if( (fbuf->iomode != DTF_IO_MODE_MEMORY) ||
            (!fbuf->explicit_match)){
            // || ( (rw_flag == DTF_WRITE) && (fbuf->writer_id != gl_my_comp_id))){

            fbuf = fbuf->next;
            continue;
        }
        /*Count for how many files we need to match ioreqs.*/
        fbuf_cnt++;
        /*Init master's db*/
        if(fbuf->mst_info->is_master_flag){
            /*Check there is no unfinished matching process*/
            if(fbuf->is_matching_flag == 1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Error: Matching for %s has not completed yet. Cannot do match all.", fbuf->file_path);
                assert(fbuf->is_matching_flag == 0);
            }
        }
        DTF_DBG(VERBOSE_DBG_LEVEL, "Will match for %s", fbuf->file_path);

        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if(ioreq->rw_flag==DTF_READ && fbuf->writer_id==gl_my_comp_id){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: cannot call dtf_match_ioreqs_all when there are uncompleted intra-component I/O requests ");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
            ioreq = ioreq->next;
        }
        /*Send ioreqs to master(s)*/
        send_ioreqs(fbuf, 0);

        /*set/reset flags*/
        fbuf->done_matching_flag = 0;
        fbuf->is_matching_flag = 1;
        fbuf = fbuf->next;
    }

    if(fbuf_cnt == 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "There are no files to match ioreqs.");

        return;
    } else
        DTF_DBG(VERBOSE_DBG_LEVEL, "Will match ioreqs for %d files", fbuf_cnt);

    /*Keep iterating over files and trying to progress I/O*/
    while(fbuf_cnt){
        tmp_cnt = 0;
        fbuf = gl_filebuf_list;
        while(fbuf != NULL){
            if( (fbuf->iomode != DTF_IO_MODE_MEMORY) ||
                (!fbuf->explicit_match)){
                //|| ( (rw_flag == DTF_WRITE) && (fbuf->writer_id != gl_my_comp_id))){

                fbuf = fbuf->next;
                continue;
            }
//TODO: what if a master receives intracomp read request from another writer? things can break?
            if(fbuf->done_matching_flag){
                tmp_cnt++;
                fbuf = fbuf->next;
                continue;
            }
            progress_io_matching();

            if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->is_master_flag)  ){
            do_matching(fbuf, 0);
            }

            fbuf = fbuf->next;
        }
        if(tmp_cnt == fbuf_cnt)
            break;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished matching I/O for all files");

    //reset all flags
    fbuf = gl_filebuf_list;
    while(fbuf != NULL){
        if( (fbuf->iomode != DTF_IO_MODE_MEMORY) ||(!fbuf->explicit_match) ){
            //|| ( (rw_flag == DTF_WRITE) && (fbuf->writer_id != gl_my_comp_id))){
            fbuf = fbuf->next;
            continue;
        }
        reset_fbuf_match(fbuf);
        fbuf = fbuf->next;
    }
}


int match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag)
{
    double t_start;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs for file %d, intracomp %d", fbuf->ncid, intracomp_io_flag);

    if(fbuf->mst_info->is_master_flag){
        //there should be no unfinished matches
        assert(fbuf->mst_info->nrranks_completed == 0);
        assert(fbuf->mst_info->nwranks_completed == 0);
    }
    t_start = MPI_Wtime();
    /*If a writer process doesn't have any io requests, it still has to
      wait for the master process to let it complete.
      If a reader process does not have any read requests,
      it notifies the master that it completed matching and returns.*/
    if(fbuf->wreq_cnt == 0 && fbuf->rreq_cnt == 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "dtf Warning: ps has no requests for file %d", fbuf->ncid);
        if(fbuf->reader_id == gl_my_comp_id){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Completed all rreqs. Notify master.");
            /*Notify master 0 that all my read io
            requests for this file have been completed*/
            int errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(errno);
            gl_stats.nmsg_sent++;
            return 0;
        }
    } else{
        DTF_DBG(VERBOSE_DBG_LEVEL, "Total %d rreqs and %d wreqs to match", fbuf->rreq_cnt, fbuf->wreq_cnt);
        send_ioreqs(fbuf, intracomp_io_flag);
    }
    //set
    fbuf->done_matching_flag = 0;
    fbuf->is_matching_flag = 1;
    /*Have to check for both things otherwise writers sometime hang
    (maybe they receive IO_CLOSE_FILE_TAG before READ_DONE_TAG by mistake?)*/

    int counter = 0;
    double t_prog_start, t_match_start, t_prog = 0, t_match=0;
    while(!fbuf->done_matching_flag){
            counter++;
            if(counter % 500 == 0){
               // sleep(0.5);
                t_prog_start = MPI_Wtime();
                progress_io_matching();
                counter = 0;
                t_prog += MPI_Wtime() - t_prog_start;

                if( (fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag  ){
                    t_match_start = MPI_Wtime();
                    do_matching(fbuf, intracomp_io_flag);
                    t_match += MPI_Wtime() - t_match_start;
                }
            }

    }

    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "tot prog %.3f, tot match %.3f", t_prog, t_match);

    gl_stats.accum_match_time += MPI_Wtime() - t_start;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Stat: Time for matching %.4f", MPI_Wtime() - t_start);

    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);

    reset_fbuf_match(fbuf);
    return 0;
}




/*writer->reader*/
static void send_data_wrt2rdr(void* buf, int bufsz)
{
    int ncid, var_id, rdr_rank, errno, i;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    dtf_var_t *var = NULL;
    size_t rofft = 0, sofft = 0;
    unsigned char *sbuf = NULL;
    size_t sbufsz = 0;
    int def_el_sz, intracomp_flag;
    unsigned char *rbuf = (unsigned char*)buf;
    MPI_Offset *start, *count;
    int nmsg_sent = 0;
    double t_recur_accum = 0, t_recur;
    int nblocks_written = 0;
    MPI_Offset *cur_coord = NULL, *new_count = NULL, *new_start = NULL, *tmp;
    MPI_Offset nelems, fit_nelems;
    double t_begin;

    rdr_rank = (int)(*(MPI_Offset*)rbuf);
    rofft += sizeof(MPI_Offset);
    intracomp_flag = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Sending data to rank %d (intracomp flag %d)", rdr_rank, intracomp_flag);
    ncid = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    assert(fbuf != NULL);

    t_begin = MPI_Wtime();
    assert(gl_msg_buf != NULL);
    sbuf = (unsigned char*)gl_msg_buf;
    sbufsz = (int)gl_conf.data_msg_size_limit;

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: Alloc data buf of sz %d", (int)sbufsz);

    while(rofft != (size_t)bufsz){

        if(sofft == 0){
            *(MPI_Offset*)sbuf = (MPI_Offset)ncid;
            sofft = sizeof(MPI_Offset);
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

        start = (MPI_Offset*)(rbuf+rofft);
        rofft += var->ndims*sizeof(MPI_Offset);
        count = (MPI_Offset*)(rbuf+rofft);
        rofft += var->ndims*sizeof(MPI_Offset);
        DTF_DBG(VERBOSE_DBG_LEVEL, "var id %d. Start -> count:", var->id);
        for(i = 0; i < var->ndims; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "       %lld --> %lld", start[i], count[i]);

        /*Find out if this portion of data fits in sbuf*/
        nelems = count[0];
        for(i = 1; i < var->ndims; i++)
            nelems *= count[i];

        {
            int min_subbl_ndims;
            size_t min_block_sz;
            size_t left_sbufsz;

            tmp = realloc(cur_coord, var->ndims*sizeof(MPI_Offset));  assert(tmp != NULL); cur_coord = tmp;

            /*For subblock that will fit in the buffer.*/
            tmp = realloc(new_count, var->ndims*sizeof(MPI_Offset)); assert(tmp != NULL);new_count = tmp;
            tmp = realloc(new_start, var->ndims*sizeof(MPI_Offset)); assert(tmp != NULL);new_start = tmp;
            for(i = 0; i < var->ndims; i++)
                new_start[i] = start[i];

            /*Define the biggest full subblock that we can fit
              if the buffer was empty*/
            size_t max_mem = sbufsz - /*ncid*/sizeof(MPI_Offset) - var->ndims*sizeof(MPI_Offset)*2 - sizeof(MPI_Offset);
            int max_nels;

            min_subbl_ndims = var->ndims;
            while(min_subbl_ndims > 0){
                max_nels = 1;
                for(i = var->ndims - min_subbl_ndims; i < var->ndims; i++)
                    max_nels *= count[i];
                if(max_nels*def_el_sz <= max_mem)
                    break;
                else
                    min_subbl_ndims--;
            }
            min_block_sz = max_nels*def_el_sz;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Will fit data in min blocks of sz %ld ndims %d out of %d dims", min_block_sz, min_subbl_ndims, var->ndims);
            assert(min_subbl_ndims>0);


            fit_nelems = 0;
            while(nelems > 0){

                if(sofft == 0){
                    *(MPI_Offset*)sbuf = (MPI_Offset)ncid;
                    sofft = sizeof(MPI_Offset);
                }
                /*How much data can we fit right now?*/
                left_sbufsz = sbufsz - sofft - var->ndims*sizeof(MPI_Offset)*2 - sizeof(MPI_Offset);

                if( left_sbufsz < min_block_sz)
                    fit_nelems = 0;
                else {
                    //adjust new_count
                    for(i = 0; i < var->ndims; i++)
                        if(i < var->ndims - min_subbl_ndims)
                            new_count[i] = 1;
                        else
                            new_count[i] = count[i];

                    fit_nelems = max_nels;

//                    if( (int)(left_sbufsz / max_nels*def_el_sz) > 1)
//                        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: copying only one min subblock but can actually fit more: %d", (int)(left_sbufsz / max_nels*def_el_sz));
                }

                if(fit_nelems == 0){
                    if(sofft == sizeof(MPI_Offset)){
                        /*The message is empty but we cannot even fit a minimum subblock*/
                        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: Cannot fit min contiguous subblock,\
                        min sz %d, increase DFT_DATA_MSG_SIZE_LIMIT (cur %d)\n", def_el_sz, (int)gl_conf.data_msg_size_limit);
                        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
                    } else {
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Send msg to %d, size %d at %.3f", rdr_rank,(int)sofft, MPI_Wtime() - gl_stats.walltime);
                        /*Send this message*/
                        MPI_Request req;
                        int dsz = (int)sofft;
                        double t_start_comm = MPI_Wtime();
                        if(intracomp_flag){
                            errno = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[gl_my_comp_id].comm, &req);
                            CHECK_MPI(errno);
                        } else {
                            errno = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm, &req);
                            CHECK_MPI(errno);
                        }
                        errno = MPI_Wait(&req, MPI_STATUS_IGNORE);
                        CHECK_MPI(errno);
                        gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                        gl_stats.nmatching_msg_sent++;
                        gl_stats.accum_msg_sz += sofft;
                        nmsg_sent++;
                        sofft = 0;
                    }
                } else {

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Will copy subblock (strt->cnt):");
                    for(i=0; i < var->ndims; i++)
                        DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld \t orig: \t %lld\t --> %lld)", new_start[i], new_count[i], start[i], count[i]);

                   /*Find the ioreq that has info about this block*/
                    int full_match;
                    ioreq = fbuf->ioreqs;
                    while(ioreq != NULL){
                        if( (ioreq->var_id == var_id) && (ioreq->rw_flag == DTF_WRITE)){
                            int match = 0;
                            full_match = 0;
                            for(i = 0; i < var->ndims; i++)
                                if( (new_start[i] >= ioreq->start[i]) && (new_start[i] < ioreq->start[i]+ioreq->count[i])){
                                    match++;
                                    if( (new_start[i] == ioreq->start[i]) && (new_count[i] == ioreq->count[i]) )
                                        full_match++;
                                }
                            if(match == var->ndims){
                                DTF_DBG(VERBOSE_DBG_LEVEL, "Matched ioreq with userbuf %p (strt->cnt): ", ioreq->user_buf);
                                for(i = 0; i < var->ndims; i++){
                                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld", ioreq->start[i], ioreq->count[i]);
                                    assert(new_start[i]+new_count[i] <= ioreq->start[i]+ioreq->count[i]);
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

                    {
                        int type_mismatch = 0;
                        if(var->dtype != ioreq->dtype){
                            /*Will only support conversion from MPI_DOUBLE to MPI_FLOAT*/
//                            assert(var->dtype == MPI_FLOAT);
//                            assert(ioreq->dtype == MPI_DOUBLE);
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
                        t_recur = MPI_Wtime();
                        *(MPI_Offset*)(sbuf+sofft) = (MPI_Offset)var_id;
                        sofft += sizeof(MPI_Offset);
                        memcpy(sbuf+sofft, new_start, var->ndims*sizeof(MPI_Offset));
                        sofft += var->ndims*sizeof(MPI_Offset);
                        memcpy(sbuf+sofft, new_count, var->ndims*sizeof(MPI_Offset));
                        sofft += var->ndims*sizeof(MPI_Offset);
                        t_recur_accum += MPI_Wtime() - t_recur;
                        nblocks_written++;

                        /*Copy data*/
                        for(i = 0; i < var->ndims; i++){
                            cur_coord[i] = new_start[i];
                        }
                        DTF_DBG(VERBOSE_DBG_LEVEL, "total data to getput %lld (of nelems %lld)", fit_nelems, nelems);

                        t_recur = MPI_Wtime();
                        if(full_match == var->ndims){

                            if(type_mismatch){
//                                int req_el_type;
//                                MPI_Type_size(ioreq->dtype, &req_el_type);
//                                assert(fit_nelems == (int)(ioreq->user_buf_sz / req_el_type));
//                                for(i = 0; i < (int)fit_nelems; i++){
//                                    ((float*)(sbuf+sofft))[i] = (float)(((double*)ioreq->user_buf)[i]);
//                                }
//                                sofft += fit_nelems*sizeof(float);
                                convertcpy(ioreq->dtype, var->dtype, ioreq->user_buf, (void*)(sbuf+sofft), (int)fit_nelems);
                            } else {
                                assert(ioreq->user_buf_sz == fit_nelems*def_el_sz);
                                memcpy(sbuf+sofft, ioreq->user_buf, ioreq->user_buf_sz);
                            }

                        } else{
                            /*read from user buffer to send buffer*/
                            recur_get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                               ioreq->count, new_start, new_count, 0,
                                               cur_coord, sbuf+sofft, DTF_READ, type_mismatch);
                        }

                        t_recur_accum += MPI_Wtime() - t_recur;
                        sofft += fit_nelems*def_el_sz+/*padding*/(fit_nelems*def_el_sz)%sizeof(MPI_Offset);
                        nelems -= fit_nelems;
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Left to get put %lld", nelems);
                        assert(nelems >= 0);
                    }

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

        }
        DTF_DBG(VERBOSE_DBG_LEVEL, "rofft %d (bufsz %d), Sofft %d", (int)rofft,bufsz, (int)sofft);
        //assert(sofft == sbufsz);

    } /*while(rofft != sbufsz)*/

    if(sofft > sizeof(MPI_Offset)){
        /*Send the last message*/
         DTF_DBG(VERBOSE_DBG_LEVEL, "Send (last) msg to %d size %d at %.3f", rdr_rank,(int)sofft, MPI_Wtime() - gl_stats.walltime);
        MPI_Request req;
        /* First send the size of the data, then the data itself*/
        int dsz = (int)sofft;
        double t_start_comm = MPI_Wtime();
        if(intracomp_flag){
            errno = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[gl_my_comp_id].comm, &req);
            CHECK_MPI(errno);
        } else {
            errno = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm, &req);
            CHECK_MPI(errno);
        }
        errno = MPI_Wait(&req, MPI_STATUS_IGNORE);
        CHECK_MPI(errno);
        gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
        gl_stats.nmatching_msg_sent++;
        gl_stats.accum_msg_sz += sofft;
        nmsg_sent++;
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: copy data to sbuf %.3f (recur+memcpy %.3f) (sent %d msgs) (%d blocks)",
            MPI_Wtime() - t_begin, t_recur_accum, nmsg_sent, nblocks_written);

    free(new_start);
    free(new_count);
    free(cur_coord);
   // dtf_free(sbuf, sbufsz); //comment out since use gl_msg_buf
}

/*writer->reader*/
static void recv_data_rdr(void* buf, int bufsz)
{
    int ncid, var_id, i;
    int def_el_sz, req_el_sz;
    dtf_var_t *var = NULL;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    ncid = (int)(*(MPI_Offset*)(chbuf));
    offt += sizeof(MPI_Offset);
    int nblocks_read = 0;
    int type_mismatch;

    MPI_Offset *start, *count;
    void *data;

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
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

        /*Find the ioreq that has info about this block*/
        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if( (ioreq->var_id == var_id) && (ioreq->rw_flag == DTF_READ)){
                int match = 0;
                for(i = 0; i < var->ndims; i++)
                    if( (start[i] >= ioreq->start[i]) && (start[i] < ioreq->start[i]+ioreq->count[i]))
                        match++;
                if(match == var->ndims){
                    for(i = 0; i < var->ndims; i++)
                        assert(start[i]+count[i]<=ioreq->start[i]+ioreq->count[i]);
                    break;
                }

            }
            ioreq = ioreq->next;
        }
        assert(ioreq != NULL);

        /*Copy data: no support for type conversion now*/
        MPI_Type_size(ioreq->dtype, &req_el_sz);
        MPI_Type_size(var->dtype, &def_el_sz);

        type_mismatch = 0;
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

        {   /*Copy data*/
            int nelems;
            int full_match = 0;
            MPI_Offset *cur_coord = (MPI_Offset*)dtf_malloc(var->ndims*sizeof(MPI_Offset));

            nelems = 1;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Getput data:");
            for(i = 0; i < var->ndims; i++){
                nelems *= count[i];
                cur_coord[i] = start[i];
                if( (ioreq->start[i] == start[i]) && (ioreq->count[i] == count[i]))
                    full_match++;
            }
//            if(req_el_sz != def_el_sz){
//                DTF_DBG(VERBOSE_DBG_LEVEL, "Dtype mismatch: def %d, req %d. Total getput datasz %d", def_el_sz, req_el_sz, def_el_sz*nelems);
//            }
           DTF_DBG(VERBOSE_DBG_LEVEL, "Will get %d elems for var %d", nelems, var->id);
           if(full_match == var->ndims){
                if(type_mismatch){

                    convertcpy(var->dtype, ioreq->dtype, data, ioreq->user_buf, nelems);
                } else {
                    assert(ioreq->user_buf_sz == nelems*def_el_sz);
                    memcpy(ioreq->user_buf, data, ioreq->user_buf_sz);
                }

            } else{
                /*Copy from rbuf to user buffer*/
                recur_get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                   ioreq->count, start, count, 0, cur_coord, data, DTF_WRITE, type_mismatch);
            }
            dtf_free(cur_coord, var->ndims*sizeof(MPI_Offset));

            offt += nelems*def_el_sz +(nelems*def_el_sz)%sizeof(MPI_Offset);
            ioreq->get_sz += (MPI_Offset)(nelems*req_el_sz);
        }

        DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, Got %d (expect %d)", ioreq->id, ioreq->var_id, (int)ioreq->get_sz, (int)ioreq->user_buf_sz);
        assert(ioreq->get_sz<=ioreq->user_buf_sz);
        if(ioreq->get_sz == ioreq->user_buf_sz){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Complete req %d (left %d)", ioreq->id, fbuf->rreq_cnt-1);

            //delete this ioreq
            delete_ioreq(fbuf, &ioreq);

            if(fbuf->rreq_cnt == 0){
                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE:Completed all rreqs");
                /*Notify my master writer rank that all my read io
                requests for this file have been completed*/
                int errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
                CHECK_MPI(errno);
                gl_stats.nmsg_sent++;

                /*Process of the reader component can set this flag now.
                  Writer component will set this flag only when the master
                  rank says so.*/
                if(fbuf->reader_id == gl_my_comp_id)
                    fbuf->done_matching_flag = 1;
            }
        }
    } //while(bufsz)

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: time to extract the data: %.3f (%d blocks)", MPI_Wtime() - t_begin, nblocks_read);
}

/*Send the file header and info about vars to the reader when the writer finishes the def mode*/
void send_file_info(file_buffer_t *fbuf, int reader_root)
{
    void *sbuf;
    MPI_Offset sbuf_sz, errno;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Will send file info to reader %d", reader_root);
    pack_file_info(fbuf, &sbuf_sz, &sbuf);
    assert(sbuf_sz > 0);
    double t_start = MPI_Wtime();
    MPI_Request req;
    errno = MPI_Isend(sbuf, (int)sbuf_sz, MPI_BYTE, reader_root, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm, &req);
    CHECK_MPI(errno);
    gl_stats.accum_comm_time += MPI_Wtime() - t_start;
    gl_stats.nmsg_sent++;
    errno = MPI_Wait(&req, MPI_STATUS_IGNORE);
    CHECK_MPI(errno);
    dtf_free(sbuf, sbuf_sz);
}

/*Function executed by master writers to notify other writers
  about something defined by mpitag (match completion or file close)*/
static void notify_workgroup(file_buffer_t *fbuf, int msgtag)
{
    MPI_Group glob_group, file_group;
    int nranks;
    int *ranks, *glob_ranks;
    int i, errno;
    int rank;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Mst %d will notify workgroup (msgtag %d) for %s", gl_my_rank, msgtag, fbuf->file_path);
    /*First, translate the ranks in the communicator which
    was used to open the file to global ranks*/
    MPI_Comm_group(gl_comps[gl_my_comp_id].comm, &glob_group);
    MPI_Comm_group(fbuf->comm, &file_group);

    nranks = fbuf->mst_info->my_workgroup_sz - 1;
    ranks = dtf_malloc(nranks*sizeof(int));
    assert(ranks != NULL);
    glob_ranks = dtf_malloc(nranks*sizeof(int));
    assert(glob_ranks != NULL);

    MPI_Comm_rank(fbuf->comm, &rank);

    for(i = 0; i < nranks; i++)
        ranks[i] = rank + i + 1;

    errno = MPI_Group_translate_ranks(file_group, nranks, ranks, glob_group, glob_ranks);
    CHECK_MPI(errno);

    for(i = 0; i < nranks; i++){
        errno = MPI_Send(&fbuf->ncid, 1, MPI_INT, glob_ranks[i], msgtag, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(errno);
        gl_stats.nmsg_sent++;
    }

    MPI_Group_free(&file_group);
    MPI_Group_free(&glob_group);
    dtf_free(ranks, nranks*sizeof(int));
    dtf_free(glob_ranks, nranks*sizeof(int));
}

void progress_io_matching()
{
    MPI_Status status;
    int comp, flag, src, errno;
    file_buffer_t *fbuf;

    int bufsz;
    void *rbuf;
    char filename[MAX_FILE_NAME];
    int i, ncid;
    double t_start = MPI_Wtime();
    double t_start_comm;
    gl_stats.nprogress_call++;
    /*first check if there are any requests for file info
      that we can process. */
    process_file_info_req_queue();

    for(comp = 0; comp < gl_ncomp; comp++){
        if( gl_comps[comp].comm == MPI_COMM_NULL){
            continue;
        }

        while(1){
            errno = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[comp].comm, &flag, &status);
            CHECK_MPI(errno);

            if(!flag){
                break;
            }
            src = status.MPI_SOURCE;

            switch(status.MPI_TAG){
                case FILE_INFO_REQ_TAG:
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv root req from rank %d (comp %s)", src, gl_comps[comp].name);
                    rbuf = dtf_malloc(MAX_FILE_NAME+2*sizeof(int));
                    assert(rbuf != NULL);
                    errno = MPI_Recv(rbuf, (int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, src, FILE_INFO_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    memcpy(filename, rbuf, MAX_FILE_NAME);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
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

                        //enqueue
                        if(gl_finfo_req_q == NULL)
                            gl_finfo_req_q = req;
                        else {
                            req->next = gl_finfo_req_q;
                            gl_finfo_req_q->prev = req;
                            gl_finfo_req_q = req;
                        }
                        break;
                    }
                    if(gl_my_rank == fbuf->root_writer){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");

                        if(fbuf->iomode == DTF_IO_MODE_FILE){
                            memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
                            /*just notify the reader that I am the root writer*/

                            errno = MPI_Send(&fbuf->ncid, 1, MPI_INT, fbuf->root_reader, FILE_INFO_TAG, gl_comps[comp].comm);
                            CHECK_MPI(errno);
                        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
                            memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
                            memcpy(&(fbuf->mst_info->nrranks_opened), (unsigned char*)rbuf+MAX_FILE_NAME+sizeof(int), sizeof(int));
                            send_file_info(fbuf, fbuf->root_reader);
                        }
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Root reader for file %s is %d", fbuf->file_path, fbuf->root_reader);

                    } else {
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Forward the request to rank %d", fbuf->root_writer);
                        /*Forward this request to the root master*/
                        double t_start = MPI_Wtime();
                        errno = MPI_Send(rbuf,(int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm);
                        CHECK_MPI(errno);

                        gl_stats.accum_comm_time += MPI_Wtime() - t_start;
                        gl_stats.nmatching_msg_sent++;
                        gl_stats.accum_msg_sz += (size_t)(MAX_FILE_NAME+2*sizeof(int));
                    }
                    dtf_free(rbuf, MAX_FILE_NAME+2*sizeof(int));
                    break;
                case IO_REQS_TAG:
                    t_start_comm = MPI_Wtime();
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Received reqs from %d, comp %d (my comp %d), bufsz %d", src, comp,gl_my_comp_id, (int)bufsz);
                    rbuf = dtf_malloc(bufsz);
                    assert(rbuf != NULL);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_REQS_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    parse_ioreqs(rbuf, bufsz, src, gl_comps[comp].comm);
                    dtf_free(rbuf, bufsz);
                    break;
                case IO_DATA_REQ_TAG:
                    t_start_comm = MPI_Wtime();
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    rbuf = dtf_malloc(bufsz);
                    assert(rbuf != NULL);

                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv data req from mst %d", src);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    send_data_wrt2rdr(rbuf, bufsz);
                    dtf_free(rbuf, bufsz);
                    break;

                case IO_DATA_TAG:
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recved data from %d", src);
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    t_start_comm = MPI_Wtime();
                    assert(bufsz>0);
                    assert(gl_msg_buf != NULL);
                    rbuf = gl_msg_buf; //dtf_malloc(bufsz);
                    assert(rbuf != NULL);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    gl_stats.accum_comm_time += MPI_Wtime() - t_start_comm;
                    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: recv data from %d", src);
                    recv_data_rdr(rbuf, bufsz);
                   // dtf_free(rbuf, bufsz);
                    break;
                case READ_DONE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, READ_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    assert(fbuf->root_writer == gl_my_rank);

                    if(comp == gl_my_comp_id){
                        fbuf->mst_info->nwranks_completed++;
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Recv read_done from writer %d (tot %d)", src, fbuf->mst_info->nwranks_completed);
                        /*shouldn't have matching going on for readers and writers at the same time*/
                        assert(fbuf->mst_info->nrranks_completed == 0);
                    } else {
                        assert(comp == fbuf->reader_id);
                        fbuf->mst_info->nrranks_completed++;
                        DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: Recv read_done from reader %d (tot %d) at %.3f",
                                src, fbuf->mst_info->nrranks_completed, MPI_Wtime() - gl_stats.walltime);
                        assert(fbuf->mst_info->nwranks_completed == 0);
                    }

                    if( ((fbuf->mst_info->nrranks_opened > 0) && (fbuf->mst_info->nrranks_completed == fbuf->mst_info->nrranks_opened)) ||
                        ((fbuf->mst_info->nwranks_opened > 0) && (fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened)) ){

                            DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: Notify masters that matching for %s is finished",fbuf->file_path);
                            for(i = 0; i < fbuf->mst_info->nmasters; i++){
                                    if(fbuf->mst_info->masters[i] == gl_my_rank)
                                        continue;
                                    DTF_DBG(VERBOSE_DBG_LEVEL, "Notify mst %d", fbuf->mst_info->masters[i]);
                                    errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], MATCH_DONE_TAG, gl_comps[gl_my_comp_id].comm);
                                    CHECK_MPI(errno);
                                    gl_stats.nmsg_sent++;
                            }

                            if(fbuf->mst_info->is_master_flag){
                                DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                                /*Tell other writer ranks that they can complete matching*/
                                notify_workgroup(fbuf, MATCH_DONE_TAG);
                                DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                                fbuf->done_matching_flag = 1;
                                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE, Done matching at %.3f", MPI_Wtime()-gl_stats.walltime);
//                                if(fbuf->rdr_closed_flag){
//                                    //DTF_DBG(VERBOSE_DBG_LEVEL, "Cleaning up all reqs and dbs");
//                                    assert(fbuf->rreq_cnt == 0);
////                                    /*Complete my own write requests
////                                      NOTE: do not delete reqs for SCALE-LETKF*/
////                                    //delete_ioreqs(fbuf);
////
////                                    if(fbuf->mst_info->is_master_flag){
////                                        /*Check that I don't have any read reqs incompleted*/
////                                        assert(fbuf->mst_info->iodb->nritems == 0);
////                                        /*Clean my iodb*/
////                                        clean_iodb(fbuf->mst_info->iodb);
////                                    }
//                                }
                            }
                    }
                    break;
                case MATCH_DONE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, MATCH_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    if(fbuf->mst_info->is_master_flag){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                        /*Tell other writer ranks that they can complete matching*/
                        notify_workgroup(fbuf, MATCH_DONE_TAG);
                    }
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                    fbuf->done_matching_flag = 1;
//                    if(fbuf->rdr_closed_flag){
////                        DTF_DBG(VERBOSE_DBG_LEVEL, "Cleaning up all reqs and dbs");
//                        assert(fbuf->rreq_cnt == 0);
//                        /*Complete my own write requests*/
//                        delete_ioreqs(fbuf);
//
//                        if(fbuf->mst_info->is_master_flag){
//                            /*Check that I don't have any read reqs incompleted*/
//                            assert(fbuf->mst_info->iodb->nritems == 0);
//                            /*Clean my iodb*/
//                            clean_iodb(fbuf->mst_info->iodb);
//                        }
                    //}
                    break;
                case IO_CLOSE_FILE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, IO_CLOSE_FILE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv close file tag for %s from %d", fbuf->file_path, src);
                    //assert(!fbuf->rdr_closed_flag); //check that haven't received that before

                    if(gl_my_rank == fbuf->root_writer){
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify other masters that readers are closing the file");
                        for(i = 0; i < fbuf->mst_info->nmasters; i++){
                            if(fbuf->mst_info->masters[i] == gl_my_rank)
                                continue;
                            errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], IO_CLOSE_FILE_TAG, gl_comps[gl_my_comp_id].comm);
                            CHECK_MPI(errno);
                            gl_stats.nmsg_sent++;
                        }
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
                    errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, ROOT_MST_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive ROOT_MST_TAG notif for %s", filename);

                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    assert(fbuf->root_writer == -1);
                    fbuf->root_writer = src;
                    break;
                case FILE_READY_TAG:
                    errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);

                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    fbuf->is_ready = 1;
                    break;
                case DONE_MULTIPLE_FLAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, DONE_MULTIPLE_FLAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
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
                default:
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: unknown tag %d", status.MPI_TAG);
                    assert(0);
            }
        }
    }

    gl_stats.t_progress += MPI_Wtime() - t_start;
}


/*function called by the writer processes*/
int init_req_match_masters(MPI_Comm comm, master_info_t *mst_info)
{

    int wg, nranks, myrank, i, errno;
    char* s = getenv("MAX_WORKGROUP_SIZE");
    MPI_Group glob_group, file_group;
    int my_master, my_master_glob;
    int *masters;

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
        mst_info->my_workgroup_sz = nranks;
        mst_info->nmasters = 1;
    } else {
        my_master = (int)(myrank/wg) * wg;
        mst_info->my_workgroup_sz = wg;
        mst_info->nmasters = (int)(nranks/wg);
        if(nranks % wg > 0){
            mst_info->nmasters++;
            if(myrank >= (mst_info->nmasters-1)*wg)
                mst_info->my_workgroup_sz = nranks % wg;
        }
    }

    if(myrank == 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Nmasters %d", mst_info->nmasters);

    mst_info->masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    assert(mst_info->masters != NULL);
    masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    assert(masters != NULL);

    masters[0] = 0;
    for(i = 1; i < mst_info->nmasters; i++){
        masters[i] = masters[i-1] + wg;
    }
    /*Translate rank from subcommunicator to the global rank in mpi_comm_world*/
    errno = MPI_Group_translate_ranks(file_group, mst_info->nmasters, masters, glob_group, mst_info->masters);
    CHECK_MPI(errno);

    /*Translate the rank of my master*/
    errno = MPI_Group_translate_ranks(file_group, 1, &my_master, glob_group, &my_master_glob);
    CHECK_MPI(errno);
    mst_info->is_master_flag = (gl_my_rank == my_master_glob) ? 1 : 0;

    if(myrank == 0){
        for(i = 0; i < mst_info->nmasters; i++)
            DTF_DBG(VERBOSE_DBG_LEVEL, "Rank %d is a master", mst_info->masters[i]);
    }
    if(mst_info->is_master_flag)
        DTF_DBG(VERBOSE_DBG_LEVEL, "My wg size %d", mst_info->my_workgroup_sz);
    dtf_free(masters, mst_info->nmasters * sizeof(int));
    MPI_Group_free(&glob_group);
    MPI_Group_free(&file_group);

    mst_info->nrranks_completed = 0;
    mst_info->nwranks_completed = 0;
    return 0;
}

