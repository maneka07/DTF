#include "pfarb_req_match.h"
#include "pfarb_common.h"
#include "pfarb_util.h"
#include "pfarb_mem.h"
#include "pfarb.h"


//static void print_ritems(read_db_item_t *items)
//{
//    read_db_item_t *ritem = items;
//    while(ritem != NULL){
//        FARB_DBG(VERBOSE_ALL_LEVEL, "Ritem for rank %d", ritem->rank);
//        ritem = ritem->next;
//    }
//}

file_info_req_q_t *gl_finfo_req_q = NULL;

/*API for handling rb_tree in write_db_item*/
void chunk_destroy(void* chunk)
{
  farb_free((write_chunk_rec_t *)chunk, sizeof(write_chunk_rec_t));
}

int chunk_cmp(const void *a, const void *b)
{
  if( ((write_chunk_rec_t*)a)->offset > ((write_chunk_rec_t*)b)->offset ) return 1;
  if( ((write_chunk_rec_t*)a)->offset + ((write_chunk_rec_t*)a)->data_sz -1 < ((write_chunk_rec_t*)b)->offset ) return -1;
  return 0;
}

void chunk_print(const void *chunk)
{
  printf("(ps %d,offt %d,sz %d)", ((write_chunk_rec_t*)chunk)->rank, (int)( ((write_chunk_rec_t*)chunk)->offset), (int)(((write_chunk_rec_t*)chunk)->data_sz));
}

void info_print(void *chunk)
{;}

void info_destroy(void *chunk)
{;}

//====================================================

void init_iodb(file_buffer_t *fbuf)
{
    fbuf->mst_info->iodb = farb_malloc(sizeof(struct ioreq_db));
    assert(fbuf->mst_info->iodb != NULL);
    fbuf->mst_info->iodb->ritems = NULL;
    fbuf->mst_info->iodb->witems = NULL;
    fbuf->mst_info->iodb->nritems = 0;
    fbuf->mst_info->iodb->updated_flag = 0;
}

static void print_read_dbitem(read_db_item_t *dbitem)
{
    read_chunk_rec_t *tmp;
    FARB_DBG(VERBOSE_ALL_LEVEL, "Read dbitem for rank %d. %d chunks:", dbitem->rank, (int)dbitem->nchunks);
    tmp = dbitem->chunks;
    while(tmp != NULL){
        FARB_DBG(VERBOSE_ALL_LEVEL, "       (var %d, %d, %d)", tmp->var_id, (int)tmp->offset, (int)tmp->data_sz);
        tmp = tmp->next;
    }
}

//static void print_write_dbitem(write_db_item_t *dbitem)
//{
//    write_chunk_rec_t *tmp;
//    FARB_DBG(VERBOSE_ALL_LEVEL, "Write dbitem for var %d.", dbitem->var_id);
//    tmp = dbitem->chunks;
//    while(tmp != NULL){
//        FARB_DBG(VERBOSE_DBG_LEVEL, "       (ps %d, off %llu, sz %llu)", tmp->rank, tmp->offset, tmp->data_sz);
//        tmp = tmp->next;
//    }
//}

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
    FARB_DBG(VERBOSE_DBG_LEVEL, "Delete io requests for file %s", fbuf->file_path);

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        delete_ioreq(fbuf, &ioreq);
        ioreq = fbuf->ioreqs;
    }
    assert(fbuf->rreq_cnt == 0);
    assert(fbuf->wreq_cnt == 0);
}

/*TODO (#9#) if file transfer type is FILE need just to send a msg FILE-READY!!!*/
static void pack_file_info(file_buffer_t *fbuf, MPI_Offset *bufsz, void **buf)
{
    farb_var_t *var;
    MPI_Offset sz = 0, offt = 0;
    unsigned char *chbuf;


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
    var = fbuf->vars;
    while(var != NULL){
        /*sz += sizeof(var->id) + sizeof(var->el_sz) +
                sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims;*/
        sz += sizeof(MPI_Offset)*3+ sizeof(MPI_Offset)*var->ndims;
        var = var->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Packing info: sz %lld", sz);
    *buf = farb_malloc(sz);
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
      FARB_DBG(VERBOSE_DBG_LEVEL, "pack %d masters", fbuf->mst_info->nmasters);
    /*list of masters*/
    memcpy(chbuf+offt, fbuf->mst_info->masters, fbuf->mst_info->nmasters+sizeof(int));
    offt += fbuf->mst_info->nmasters*sizeof(MPI_Offset); //sizeof(int) + padding for MPI_Offset
    /*number of vars*/
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->var_cnt;
    offt += sizeof(MPI_Offset);
    FARB_DBG(VERBOSE_DBG_LEVEL, "pack %d vars", fbuf->var_cnt);
    /*vars*/
    var = fbuf->vars;
    while(var != NULL){
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->id;
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)mpitype2int(var->dtype);
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->ndims;
        offt += sizeof(MPI_Offset);
        memcpy((void*)(chbuf+offt), (void*)var->shape, sizeof(MPI_Offset)*var->ndims);
        offt += sizeof(MPI_Offset)*var->ndims;
        var = var->next;
    }
    FARB_DBG(VERBOSE_ALL_LEVEL, "offt %lld", offt);
    assert(offt == sz);
    *bufsz = sz;
}

void unpack_file_info(MPI_Offset bufsz, void *buf)
{
    int i, varid, var_cnt;
    file_buffer_t *fbuf;
    farb_var_t *var;
    MPI_Offset offt = 0;
    int type;
    unsigned char *chbuf = (unsigned char*)buf;
    char filename[MAX_FILE_NAME];
    FARB_DBG(VERBOSE_DBG_LEVEL, "Start unpackinf file info for of sz %d", (int)bufsz);
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
    FARB_DBG(VERBOSE_DBG_LEVEL, "unpack filename %s", filename);
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);
    /*ncid*/
    fbuf->ncid = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    /*header size*/
    fbuf->hdr_sz = *((MPI_Offset*)(chbuf+offt));
    offt += sizeof(MPI_Offset);
    /*header*/
    fbuf->header = farb_malloc(fbuf->hdr_sz);
    assert(fbuf->header != NULL);
    memcpy(fbuf->header, chbuf+offt, fbuf->hdr_sz);
    offt += fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset);
    /*number of masters*/
    fbuf->mst_info->nmasters = (int)(*((MPI_Offset*)(chbuf+offt)));
    FARB_DBG(VERBOSE_DBG_LEVEL, "unpack %d masters", fbuf->mst_info->nmasters);
    offt += sizeof(MPI_Offset);
    /*list of masters*/
    fbuf->mst_info->masters = farb_malloc(fbuf->mst_info->nmasters*sizeof(int));
    assert(fbuf->mst_info->masters != NULL);
    memcpy(fbuf->mst_info->masters, chbuf+offt, fbuf->mst_info->nmasters*sizeof(int));
    offt += fbuf->mst_info->nmasters*sizeof(MPI_Offset);
    //init root master
    fbuf->root_writer = fbuf->mst_info->masters[0];
    /*number of vars*/
    var_cnt = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    FARB_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", var_cnt);
    /*vars*/
    for(i = 0; i < var_cnt; i++){
        varid = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        var = find_var(fbuf->vars, varid);
        assert(var == NULL);
        var = new_var(varid, 0, 0, NULL);
        type = (int)(*((MPI_Offset*)(chbuf+offt)));
        var->dtype = int2mpitype(type);
        offt+=sizeof(MPI_Offset);
        var->ndims = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);

        if(var->ndims > 0){
            var->shape = farb_malloc(var->ndims*sizeof(MPI_Offset));
            assert(var->shape != NULL);
            memcpy((void*)var->shape, chbuf+offt, sizeof(MPI_Offset)*var->ndims);
            offt += sizeof(MPI_Offset)*var->ndims;
        } else
            var->shape = NULL;

        add_var(&(fbuf->vars), var);
        fbuf->var_cnt++;
    }
    assert(offt == bufsz);
    FARB_DBG(VERBOSE_ALL_LEVEL, "Finished unpacking");
}

/*Master rank tries to match as many memchunks in the read request as it can
and forwards corresponding part of the read request to those writer ranks that
hold the data */
static void do_matching(file_buffer_t *fbuf, int intracomp_io_flag)
{
    int mlc_ranks = 4;
    int mlc_buf   = 512;
    int i;
    read_chunk_rec_t  *rchunk;
    write_chunk_rec_t  match_chunk;
    int rank_idx, var_id;
    int nwriters, allocd_nwriters;
    write_db_item_t *witem = NULL;
    read_db_item_t *ritem;
    int *writers;
    unsigned char **sbuf;
    int *bufsz;
    size_t *offt;
    int matched_rank;
    MPI_Offset matched_offt, matched_sz;


    double t_start, t_send = 0;

    //if(!fbuf->mst_info->iodb->updated_flag) //no new info since last time matching was done, ignore
      //  return;
    fbuf->mst_info->iodb->updated_flag = 0; //reset

    t_start = MPI_Wtime();

    writers = (int*)farb_malloc(mlc_ranks*sizeof(int));
    assert(writers != NULL);
    sbuf = (unsigned char**)farb_malloc(mlc_ranks*sizeof(unsigned char*));
    assert(sbuf != NULL);
    bufsz = (int*)farb_malloc(mlc_ranks*sizeof(int));
    assert(bufsz != NULL);
    offt = (size_t*)farb_malloc(mlc_ranks*sizeof(size_t));
    assert(offt != NULL);

    allocd_nwriters = mlc_ranks;

    for(i = 0; i < mlc_ranks; i++){
        sbuf[i] = NULL;
        bufsz[i] = 0;
        offt[i] = 0;
    }

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
            FARB_DBG(VERBOSE_DBG_LEVEL, "Writer has unfinished intracomp rreqs for rank %d, chunks are:", ritem->rank);
            print_read_dbitem(ritem);
            assert(0);
            ritem = ritem->next;
            continue;
        }


        nwriters = 0;

        for(i = 0; i < allocd_nwriters; i++){
            //sbuf[i] = NULL;
            //bufsz[i] = 0;
            offt[i] = 0;
            //writers[i] = -1;
        }
        if(ritem->comm == gl_comps[gl_my_comp_id].comm)
            FARB_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in my comp, chunks before matching:", ritem->rank);
        else
            FARB_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in other comp, chunks before matching:", ritem->rank);
        print_read_dbitem(ritem);

        rchunk = ritem->chunks;
        while(rchunk != NULL){

            var_id = rchunk->var_id;

            if( (witem == NULL) || (witem->var_id != var_id)){
                //find the write record for this var_id
                witem = fbuf->mst_info->iodb->witems;
                while(witem != NULL){
                    if(witem->var_id == var_id)
                        break;
                    witem = witem->next;
                }
                if(witem == NULL){
                    rchunk = rchunk->next;
                    /*No match right now*/
                    continue;
                }
            }
            FARB_DBG(VERBOSE_ALL_LEVEL, "write record for this var:");
//            print_write_dbitem(witem);

            rb_red_blk_node *wchunk;
            while(rchunk->data_sz){
                match_chunk.offset = rchunk->offset;
                match_chunk.data_sz = rchunk->data_sz;
                wchunk = RBExactQuery(witem->chunks, &match_chunk);
                if(wchunk == NULL){
                    //couldn't match
                    //rchunk = rchunk->next;
                    break;
                }

                matched_rank = ((write_chunk_rec_t*)wchunk->key)->rank;
                matched_offt = match_chunk.offset;
                if(match_chunk.offset+match_chunk.data_sz >
                    ((write_chunk_rec_t*)wchunk->key)->offset + ((write_chunk_rec_t*)wchunk->key)->data_sz){
                    //matched only part of the chunk
                    matched_sz = ((write_chunk_rec_t*)wchunk->key)->offset +
                                 ((write_chunk_rec_t*)wchunk->key)->data_sz -
                                 match_chunk.offset;
                } else {
                    matched_sz = match_chunk.data_sz;
                }

                {
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
                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            writers = (int*)tmp;

                            tmp = realloc((void*)offt, (nwriters+mlc_ranks)*sizeof(size_t));
                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            offt = (size_t*)tmp;

                            tmp = realloc((void*)bufsz, (nwriters+mlc_ranks)*sizeof(int));
                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
                            assert(tmp != NULL);
                            bufsz = (int*)tmp;

                            tmp1 = realloc(sbuf, (nwriters+mlc_ranks)*sizeof(unsigned char*));
                            gl_conf.malloc_size += mlc_ranks*sizeof(unsigned char*);
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
                            gl_conf.malloc_size += mlc_buf;
                        }

                        *(MPI_Offset*)(sbuf[rank_idx]) = (MPI_Offset)ritem->rank; //to whom writer should send the data
                        offt[rank_idx] += sizeof(MPI_Offset);
                        if(intracomp_io_flag)
                            FARB_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to another w");
                        else
                           FARB_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to reader");
                        *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)intracomp_io_flag; //intra- or inter- component?
                        offt[rank_idx] += sizeof(MPI_Offset);
                        *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)fbuf->ncid;  //data from what file
                        offt[rank_idx] += sizeof(MPI_Offset);
                     }
                     /*Save infor about the mem chunk*/
                    //FARB_DBG(VERBOSE_ALL_LEVEL, "rank idx %d, nranks %d", rank_idx, nwriters);
                    if(offt[rank_idx] + sizeof(MPI_Offset)*3 > bufsz[rank_idx]){
                        unsigned char *tmp;
                        tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
                        assert(tmp != NULL);
                        sbuf[rank_idx] = tmp;
                        bufsz[rank_idx] += mlc_buf;
                        gl_conf.malloc_size += mlc_buf;
                    }

                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)rchunk->var_id;
                    offt[rank_idx] += sizeof(MPI_Offset);
                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = matched_offt;
                    offt[rank_idx] += sizeof(MPI_Offset);
                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = matched_sz;
                    offt[rank_idx] += sizeof(MPI_Offset);

                    FARB_DBG(VERBOSE_ALL_LEVEL, "Will ask %d for data (%d, %d)", matched_rank, (int)matched_offt, (int)matched_sz);
                }

                rchunk->offset += matched_offt;
                rchunk->data_sz -= matched_sz;
            }

            if(rchunk->data_sz == 0){
                //matched all
                FARB_DBG(VERBOSE_ALL_LEVEL, "Matched all chunk (%d, %d). Will delete (nchunk %d)", (int)rchunk->offset, (int)rchunk->data_sz, (int)ritem->nchunks);
                /*delete this chunk*/
                read_chunk_rec_t *tmp = rchunk;
                if(rchunk == ritem->last)
                    ritem->last = ritem->last->prev;
                if(rchunk == ritem->chunks)
                    ritem->chunks = ritem->chunks->next;
                if(rchunk->next != NULL)
                    rchunk->next->prev = rchunk->prev;
                if(rchunk->prev != NULL)
                    rchunk->prev->next = rchunk->next;

                rchunk = rchunk->next;
                farb_free(tmp, sizeof(read_chunk_rec_t));
                ritem->nchunks--;
                continue;
            }

            rchunk = rchunk->next;
        }
            {
            //TODO any better way than comparing each rchunk with each wchunk every time?
//            wchunk = witem->chunks;
//            while(wchunk != NULL){
//                matchsz = 0;
//                if( (rchunk->offset >= wchunk->offset) && (rchunk->offset < wchunk->offset + wchunk->data_sz)){
//                    start_offt = rchunk->offset;
//                    if( start_offt + rchunk->data_sz <= wchunk->offset + wchunk->data_sz)
//                        matchsz = rchunk->data_sz;
//                    else
//                        matchsz = wchunk->offset + wchunk->data_sz - start_offt;
//                } else if( (wchunk->offset > rchunk->offset) && (wchunk->offset < rchunk->offset + rchunk->data_sz)){
//                    start_offt = wchunk->offset;
//
//                    if(wchunk->offset + wchunk->data_sz >= rchunk->offset + rchunk->data_sz)
//                        matchsz = rchunk->offset + rchunk->data_sz - start_offt;
//                    else
//                        matchsz = wchunk->data_sz;
//                } else {
//                    wchunk = wchunk->next;
//                    continue;
//                }
//
//                /*Find send buffer for this rank*/
//                if(nwriters == 0){
//                    writers[nwriters] = wchunk->rank;
//                    rank_idx = 0;
//                    nwriters++;
//                } else {
//                    for(i = 0; i < nwriters; i++)
//                        if(writers[i] == wchunk->rank)
//                            break;
//                    if(i == nwriters){
//                        /*add new send buffer*/
//                        if(nwriters % mlc_ranks == 0 ){
//                            void *tmp;
//                            unsigned char **tmp1;
//                            tmp = realloc((void*)writers, (nwriters+mlc_ranks)*sizeof(int));
//                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
//                            assert(tmp != NULL);
//                            writers = (int*)tmp;
//
//                            tmp = realloc((void*)offt, (nwriters+mlc_ranks)*sizeof(size_t));
//                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
//                            assert(tmp != NULL);
//                            offt = (size_t*)tmp;
//
//                            tmp = realloc((void*)bufsz, (nwriters+mlc_ranks)*sizeof(int));
//                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
//                            assert(tmp != NULL);
//                            bufsz = (int*)tmp;
//
//                            tmp1 = realloc(sbuf, (nwriters+mlc_ranks)*sizeof(unsigned char*));
//                            gl_conf.malloc_size += mlc_ranks*sizeof(int);
//                            assert(tmp1 != NULL);
//                            sbuf = tmp1;
//                        }
//                        writers[nwriters] = wchunk->rank;
//                        offt[nwriters] = 0;
//                        bufsz[nwriters] = 0;
//                        sbuf[nwriters] = NULL;
//                        rank_idx = nwriters;
//                        nwriters++;
//                    } else
//                        rank_idx = i;
//                }
//                if(offt[rank_idx] == 0) {
//                    unsigned char *tmp;
//                    tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
//                    gl_conf.malloc_size += mlc_buf;
//                    assert(tmp != NULL);
//                    sbuf[rank_idx] = tmp;
//                    bufsz[rank_idx] += mlc_buf;
//                    *(MPI_Offset*)(sbuf[rank_idx]) = (MPI_Offset)ritem->rank; //to whom writer should send the data
//                    offt[rank_idx] += sizeof(MPI_Offset);
//                    if(intracomp_io_flag)
//                        FARB_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to another w");
//                    else
//                       FARB_DBG(VERBOSE_DBG_LEVEL, "Tell w to send the data to reader");
//                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)intracomp_io_flag; //intra- or inter- component?
//                    offt[rank_idx] += sizeof(MPI_Offset);
//                    *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)fbuf->ncid;  //data from what file
//                    offt[rank_idx] += sizeof(MPI_Offset);
//                }
//                //FARB_DBG(VERBOSE_ALL_LEVEL, "rank idx %d, nranks %d", rank_idx, nwriters);
//                if(offt[rank_idx] + sizeof(MPI_Offset)*3 > bufsz[rank_idx]){
//                    unsigned char *tmp;
//                    tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + mlc_buf);
//                    gl_conf.malloc_size += mlc_buf;
//                    assert(tmp != NULL);
//                    sbuf[rank_idx] = tmp;
//                    bufsz[rank_idx] += mlc_buf;
//                }
//
//                *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)rchunk->var_id;
//                offt[rank_idx] += sizeof(MPI_Offset);
//                *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = start_offt;
//                offt[rank_idx] += sizeof(MPI_Offset);
//                *(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = matchsz;
//                offt[rank_idx] += sizeof(MPI_Offset);
//
//                FARB_DBG(VERBOSE_ALL_LEVEL, "Will ask %d for data (%d, %d)", wchunk->rank, (int)start_offt, (int)matchsz);
//                /*Check if all the data of this chunk has been matched*/
//                if(matchsz == rchunk->data_sz){
//                    break;
//                } else if(matchsz > 0){
//                    /* For now, we managed to match only a part of memory chunk*/
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Matched part of req: (%d, %d) of (%d, %d)", (int)start_offt, (int)matchsz, (int)rchunk->offset, (int)rchunk->data_sz);
//                    if(start_offt == rchunk->offset){
//                        rchunk->offset = start_offt + matchsz;
//                        rchunk->data_sz -= matchsz;
//                    } else if(start_offt + matchsz == rchunk->offset + rchunk->data_sz){
//                        rchunk->data_sz -= matchsz;
//                    } else {
//                        FARB_DBG(VERBOSE_DBG_LEVEL, "Matched in the middle of a chunk");
//                        /*We matched data in the middle of the chunk,
//                        hence, we have to split it into 2 chunks now*/
//                        read_chunk_rec_t *chunk = farb_malloc(sizeof(read_chunk_rec_t));
//                        assert(chunk != NULL);
//                        chunk->var_id = rchunk->var_id;
//                        chunk->offset = start_offt + matchsz;
//                        chunk->data_sz = rchunk->offset + rchunk->data_sz - chunk->offset;
//                        rchunk->data_sz -= matchsz + chunk->data_sz;
//                        /*Instert the new chunk and skip it*/
//                        chunk->next = rchunk->next;
//                        if(chunk->next != NULL)
//                            chunk->next->prev = chunk;
//                        rchunk->next = chunk;
//                        chunk->prev = rchunk;
//                        ritem->nchunks++;
//                    }
//
//                    print_read_dbitem(ritem);
//                }
//                wchunk = wchunk->next;
//            }
//
//            if(matchsz == rchunk->data_sz){
//                FARB_DBG(VERBOSE_ALL_LEVEL, "Matched all chunk (%d, %d). Will delete (nchunk %d)", (int)rchunk->offset, (int)rchunk->data_sz, (int)ritem->nchunks);
//                /*delete this chunk*/
//                read_chunk_rec_t *tmp = rchunk;
//                if(rchunk == ritem->chunks)
//                    ritem->chunks = ritem->chunks->next;
//                if(rchunk->next != NULL)
//                    rchunk->next->prev = NULL;
//                if(rchunk->prev != NULL)
//                    rchunk->prev->next = rchunk->next;
//
//                rchunk = rchunk->next;
//                farb_free(tmp, sizeof(read_chunk_rec_t));
//                ritem->nchunks--;
//                continue;
//            }
//            rchunk = rchunk->next;
//        }
            }

        if(nwriters > 0){
            int errno;
            int completed;
            MPI_Request *sreqs = (MPI_Request*)farb_malloc(nwriters*sizeof(MPI_Request));
            assert(sreqs != NULL);

            double t_start_send = MPI_Wtime();
            for(i = 0; i < nwriters; i++){
                assert(offt[i] > 0);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Send data req to wrt %d", writers[i]);
                errno = MPI_Isend((void*)sbuf[i], offt[i], MPI_BYTE, writers[i], IO_DATA_REQ_TAG, gl_comps[gl_my_comp_id].comm, &sreqs[i]);
                CHECK_MPI(errno);
                gl_stats.nmatching_msg_sent++;
                gl_stats.accum_msg_sz += offt[i];
            }

//            errno = MPI_Waitall(nwriters, sreqs, MPI_STATUSES_IGNORE);
//            CHECK_MPI(errno);
            completed = 0;
            //FARB_DBG(VERBOSE_ALL_LEVEL, "Start waiting (rreq from rank %d)", ritem->rank);
            while(!completed){
                progress_io_matching();
                errno = MPI_Testall(nwriters, sreqs, &completed, MPI_STATUSES_IGNORE);
                CHECK_MPI(errno);
            }
            t_send += MPI_Wtime() - t_start_send;
            //FARB_DBG(VERBOSE_ALL_LEVEL, "Finish waiting (rreq from rank %d)", ritem->rank);
            farb_free(sreqs, nwriters*sizeof(MPI_Request));
//            for(i = 0; i < nwriters; i++)
//                farb_free(sbuf[i], offt[i]);

            FARB_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for rank %d", ritem->rank);
            print_read_dbitem(ritem);
        } //else {
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Did not manage to match anything for rank %d", ritem->rank);
//        }

        /*If we matched all chunks for this rank, then delete this ritem*/
        if(ritem->nchunks == 0){
            read_db_item_t *tmp = ritem;
            FARB_DBG(VERBOSE_DBG_LEVEL, "Matched all. Delete ritem of rank %d (left ritems %d). ", ritem->rank,  (int)(fbuf->mst_info->iodb->nritems - 1));
            //FARB_DBG(VERBOSE_ALL_LEVEL, "Ritems before:");
            //print_ritems(fbuf->mst_info->iodb->ritems);
            if(ritem == fbuf->mst_info->iodb->ritems)
                fbuf->mst_info->iodb->ritems = ritem->next;//fbuf->mst_info->iodb->ritems->next;
            if(ritem->prev != NULL)
                ritem->prev->next = ritem->next;
            if(ritem->next != NULL)
                ritem->next->prev = ritem->prev;
            ritem = ritem->next;
            farb_free(tmp, sizeof(read_db_item_t));
            fbuf->mst_info->iodb->nritems--;
        } else {
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Not all chunks for rreq from rank %d have been matched (%d left)", ritem->rank, (int)ritem->nchunks);
//            FARB_DBG(VERBOSE_ALL_LEVEL, "Chunks left:");
            print_read_dbitem(ritem);
            ritem = ritem->next;
        }
    }

    for(i = 0; i < allocd_nwriters; i++){
        if(sbuf[i] != NULL)
            farb_free(sbuf[i], (size_t)bufsz[i]);
    }
    /*dealloc stuff*/
    farb_free(writers, allocd_nwriters*sizeof(int));
    farb_free(offt, allocd_nwriters*sizeof(size_t));
    farb_free(bufsz, allocd_nwriters*sizeof(int));
    farb_free(sbuf, allocd_nwriters*sizeof(unsigned char*));

    gl_stats.ndb_match++;
    gl_stats.accum_db_match_time += MPI_Wtime() - t_start;
    FARB_DBG(VERBOSE_ERROR_LEVEL, "Stat: Time to match db reqs %.4f, time to send %.4f", MPI_Wtime() - t_start, t_send);

}

void clean_iodb(ioreq_db_t *iodb)
{
    write_db_item_t *witem;
    read_db_item_t *ritem;
    read_chunk_rec_t *rrec;

    size_t old_mem_sz;

    witem = iodb->witems;
    while(witem != NULL){
//        wrec = witem->chunks;
//        while(wrec != NULL){
//            witem->chunks = witem->chunks->next;
//            farb_free(wrec, sizeof(write_chunk_rec_t));
//            wrec = witem->chunks;
//        }

        old_mem_sz = gl_conf.malloc_size;
        RBTreeDestroy(witem->chunks);
        //check that no mem leakage
        FARB_DBG(VERBOSE_DBG_LEVEL, "mem before destroy tree %lu, after %lu, should differ by %lu",
                 old_mem_sz, gl_conf.malloc_size, (size_t)witem->nchunks*sizeof(write_chunk_rec_t));
        assert(old_mem_sz - gl_conf.malloc_size  == (size_t)witem->nchunks*sizeof(write_chunk_rec_t));
        iodb->witems = iodb->witems->next;
        farb_free(witem, sizeof(write_db_item_t));
        witem = iodb->witems;
    }

    ritem = iodb->ritems;
    while(ritem != NULL){
        rrec = ritem->chunks;
        while(rrec != NULL){
            ritem->chunks = ritem->chunks->next;
            farb_free(rrec, sizeof(read_chunk_rec_t));
            rrec = ritem->chunks;
        }

        iodb->ritems = iodb->ritems->next;
        farb_free(ritem, sizeof(read_db_item_t));
        ritem = iodb->ritems;
        iodb->nritems--;
    }

    iodb->witems = NULL;
    iodb->ritems = NULL;
    assert(iodb->nritems == 0);
}

//static write_chunk_rec_t *new_write_chunk_rec(int rank, MPI_Offset offset, MPI_Offset data_sz)
//{
//    write_chunk_rec_t *chunk = farb_malloc(sizeof(write_chunk_rec_t));
//    assert(chunk != NULL);
//    chunk->rank = rank;
//    chunk->offset = offset;
//    chunk->data_sz = data_sz;
//    chunk->next = NULL;
//    chunk->prev = NULL;
//    return chunk;
//}


//static void parse_write_ioreq(file_buffer_t *fbuf, farb_var_t *var, void *buf, int bufsz, int rank, MPI_Comm comm);
//static void parse_read_ioreq(file_buffer_t *fbuf, farb_var_t *var, void *buf, int bufsz, int rank, MPI_Comm comm);

static void parse_ioreqs(void *buf, int bufsz, int rank, MPI_Comm comm)
{
    int ncid, var_id, rw_flag, nchunks;
    file_buffer_t *fbuf;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    ncid = (int)(*((MPI_Offset*)chbuf));
    offt += sizeof(MPI_Offset);
    double t_start = MPI_Wtime();

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    assert(fbuf != NULL);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Start parsing reqs for file %s", fbuf->file_path);
    if(comm == gl_comps[fbuf->reader_id].comm)
        FARB_DBG(VERBOSE_DBG_LEVEL, "Reqs are from reader");
    else {
        assert(comm == gl_comps[fbuf->writer_id].comm);
        FARB_DBG(VERBOSE_DBG_LEVEL, "Req are from writer");
    }

    while(offt != (size_t)bufsz){
        rw_flag = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        nchunks = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);

        if(nchunks == 0)
            continue;

        if(rw_flag == FARB_READ){
            FARB_DBG(VERBOSE_DBG_LEVEL, "its a rreq");
            read_db_item_t *dbitem;
            read_chunk_rec_t *chunk;
            /*Find corresponding record in the database*/
            dbitem = fbuf->mst_info->iodb->ritems;
            while(dbitem != NULL){
                if( (dbitem->rank == rank) && (dbitem->comm == comm) )
                    break;
                dbitem = dbitem->next;
            }

            if(dbitem == NULL){
                dbitem = (read_db_item_t*)farb_malloc(sizeof(read_db_item_t));
                assert(dbitem != NULL);
                dbitem->rank = rank;
                dbitem->comm = comm;
                dbitem->next = NULL;
                dbitem->prev = NULL;
                dbitem->chunks = NULL;
                dbitem->last = NULL;
                dbitem->nchunks = nchunks;
                //enqueue
                if(fbuf->mst_info->iodb->ritems == NULL)
                    fbuf->mst_info->iodb->ritems = dbitem;
                else{
                    dbitem->next = fbuf->mst_info->iodb->ritems;
                    fbuf->mst_info->iodb->ritems->prev = dbitem;
                    fbuf->mst_info->iodb->ritems = dbitem;
                }
                fbuf->mst_info->iodb->nritems++;
            } else
                dbitem->nchunks += nchunks;
            FARB_DBG(VERBOSE_DBG_LEVEL, "ritems %d",(int)fbuf->mst_info->iodb->nritems);
            FARB_DBG(VERBOSE_DBG_LEVEL, "Add %d chunks", nchunks);
            //add to read items
            while(nchunks){

                chunk = (read_chunk_rec_t*)farb_malloc(sizeof(read_chunk_rec_t));
                assert(chunk != NULL);
                chunk->var_id = var_id;
                chunk->offset = *(MPI_Offset*)(chbuf+offt);
                offt+=sizeof(MPI_Offset);
                chunk->data_sz = *(MPI_Offset*)(chbuf+offt);
                offt+=sizeof(MPI_Offset);
                chunk->next = NULL;
                chunk->prev = NULL;

                if(dbitem->chunks == NULL){
                    dbitem->chunks = chunk;
                    dbitem->last = chunk;
                } else{
                    dbitem->last->next = chunk;
                    chunk->prev = dbitem->last;
                    dbitem->last = chunk;
                }
                nchunks--;
            }
        } else {
            FARB_DBG(VERBOSE_DBG_LEVEL, "its a wreq");
            write_chunk_rec_t *chunk;
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
                dbitem = (write_db_item_t*)farb_malloc(sizeof(write_db_item_t));
                assert(dbitem != NULL);
                dbitem->var_id = var_id;
                dbitem->next = NULL;
                dbitem->chunks = RBTreeCreate(chunk_cmp, chunk_destroy, info_destroy, chunk_print, info_print);
                assert(dbitem->chunks != NULL);
                //enqueue
                if(fbuf->mst_info->iodb->witems == NULL)
                    fbuf->mst_info->iodb->witems = dbitem;
                else{
                    dbitem->next = fbuf->mst_info->iodb->witems;
                    fbuf->mst_info->iodb->witems = dbitem;
                }
                dbitem->nchunks = (MPI_Offset)nchunks;
            } else
                dbitem->nchunks += (MPI_Offset)nchunks;

            FARB_DBG(VERBOSE_DBG_LEVEL, "Add %d chunks", nchunks);
            //add to write items
            while(nchunks){
                chunk = (write_chunk_rec_t*)farb_malloc(sizeof(write_chunk_rec_t));
                assert(chunk != NULL);
                chunk->offset = *(MPI_Offset*)(chbuf+offt);
                offt+=sizeof(MPI_Offset);
                chunk->data_sz = *(MPI_Offset*)(chbuf+offt);
                offt+=sizeof(MPI_Offset);
                chunk->rank = rank;
                //TODO the library doesnt detect if there is overlapping!!
                RBTreeInsert(dbitem->chunks, chunk, 0);
                nchunks--;
            }
        }
        FARB_DBG(VERBOSE_DBG_LEVEL, "cur offt %lu, bufsz %d", offt, bufsz);
    }
    assert(offt == (size_t)bufsz);
    fbuf->mst_info->iodb->updated_flag = 1;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Stat: time to parse reqs %.4f", MPI_Wtime()-t_start);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished parsing reqs. (mem %lu)", gl_conf.malloc_size);

}


//static void parse_ioreq(void *buf, int bufsz, int rank, MPI_Comm comm)
//{
//    int ncid, var_id, rw_flag;
//    file_buffer_t *fbuf;
//    farb_var_t *var;
//    size_t offt = 0;
//    unsigned char *chbuf = (unsigned char*)buf;
//    ncid = (int)(*((MPI_Offset*)chbuf));
//    offt += sizeof(MPI_Offset);
//    rw_flag = (int)(*(MPI_Offset*)(chbuf+offt));
//    offt += sizeof(MPI_Offset);
//    var_id = (int)(*(MPI_Offset*)(chbuf+offt));
//    offt += sizeof(MPI_Offset);
//
//    if(rw_flag == FARB_READ)
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Master recv read req from %d for ncid %d, var %d", rank, ncid, var_id);
//    else
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Master recv write req from %d for ncid %d, var %d", rank, ncid, var_id);
//
//    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
//    assert(fbuf != NULL);
//
//    if(comm == gl_comps[fbuf->reader_id].comm)
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Req is from reader");
//    else {
//        assert(comm == gl_comps[fbuf->writer_id].comm);
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Req is from writer");
//    }
//
//    var = find_var(fbuf->vars, var_id);
//    assert(var != NULL);
//
//    bufsz -= (int)offt;
//
//    if(rw_flag == FARB_READ)
//        parse_read_ioreq(fbuf, var, chbuf+offt, bufsz, rank, comm);
//    else
//        parse_write_ioreq(fbuf, var, chbuf+offt, bufsz, rank, comm);
//    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished parsing. (mem %lu)", gl_conf.malloc_size);
//
//}

/*
    update the database with write io request
*/
//static void parse_write_ioreq(file_buffer_t *fbuf, farb_var_t *var, void *buf, int bufsz, int rank, MPI_Comm comm)
//{
//    write_db_item_t *dbitem;
//    write_chunk_rec_t *new_chunk;
//
//    int mst1, mst2, mst;
//    MPI_Offset offset, data_sz, offt_range;
//    unsigned char *chbuf = (unsigned char*)buf;
//    size_t offt = 0;
//
//    /*Allow write requests only from the writer*/
//    assert(comm == gl_comps[fbuf->writer_id].comm);
//
//    /*Find corresponding record in the database*/
//    assert(fbuf->mst_info->iodb!=NULL);
//    dbitem = fbuf->mst_info->iodb->witems;
//    while(dbitem != NULL){
//        if(dbitem->var_id == var->id)
//            break;
//        dbitem = dbitem->next;
//    }
//
//    if(dbitem == NULL){
//        dbitem = (write_db_item_t*)farb_malloc(sizeof(write_db_item_t));
//        assert(dbitem != NULL);
//        dbitem->var_id = var->id;
//        dbitem->next = NULL;
//        dbitem->chunks = NULL;
//        dbitem->last = NULL;
//        dbitem->memchunks = NULL;
//
//        //enqueue
//        if(fbuf->mst_info->iodb->witems == NULL)
//            fbuf->mst_info->iodb->witems = dbitem;
//        else{
//            dbitem->next = fbuf->mst_info->iodb->witems;
//            fbuf->mst_info->iodb->witems = dbitem;
//        }
//    }
//
//    /*add mem chunks to the list for future matching.
//      Pick only those mem chunks whose offsets are within
//      the range I'm responsible for*/
//
//    if(has_unlim_dim(var))
//        offt_range = 0;
//    else {
//        int el_sz;
//        MPI_Type_size(var->dtype, &el_sz);
//        FARB_DBG(VERBOSE_ALL_LEVEL, "last idx %d, offt %d", (int)last_1d_index(var->ndims, var->shape), (int)(last_1d_index(var->ndims, var->shape)*el_sz));
//        offt_range = (MPI_Offset)((last_1d_index(var->ndims, var->shape) + 1)/fbuf->mst_info->nmasters);
//        offt_range *= el_sz;
//    }
//
//    FARB_DBG(VERBOSE_ALL_LEVEL, "offt_range %d", (int)offt_range);
//
//    FARB_DBG(VERBOSE_DBG_LEVEL, "List before update (var %d)", var->id);
//    print_write_dbitem(dbitem);
//    while(offt != (size_t)bufsz){
//
//        offset = *(MPI_Offset*)(chbuf+offt);
//        offt+=sizeof(MPI_Offset);
//        data_sz = *(MPI_Offset*)(chbuf+offt);
//        offt+=sizeof(MPI_Offset);
//
//        if(offt_range == 0){//Number of elements in the variable is less than then number of masters
//            mst1 = mst2 = 0;
//        } else {
//            mst1 = (int)(offset/offt_range);
//            if(mst1 >= fbuf->mst_info->nmasters)
//                mst1 = fbuf->mst_info->nmasters - 1;
//            mst2 = (int)( (offset+data_sz-1)/offt_range);
//            if(mst2 >= fbuf->mst_info->nmasters)
//                mst2 = fbuf->mst_info->nmasters - 1;
//        }
//
//        FARB_DBG(VERBOSE_ALL_LEVEL, "--> (%d, %d): mst1 %d mst2 %d", (int)offset, (int)data_sz, mst1, mst2);
//        if( (fbuf->mst_info->masters[mst1] != gl_my_rank) && (fbuf->mst_info->masters[mst2] != gl_my_rank))
//            /*This master is not responsible for this data*/
//            continue;
//        else if(fbuf->mst_info->masters[mst1] != gl_my_rank){
//            for(mst = 0; mst < fbuf->mst_info->nmasters; mst++)
//                if(fbuf->mst_info->masters[mst] == gl_my_rank)
//                    break;
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Changing offset from %llu to %llu, datasz from %llu to %llu",
//                     offset, offt_range*mst, data_sz, data_sz - (offt_range*mst - offset));
//            data_sz -= offt_range*mst - offset;
//            offset = offt_range*mst;
//
//        } else if (fbuf->mst_info->masters[mst2] != gl_my_rank){
//            for(mst = 0; mst < fbuf->mst_info->nmasters; mst++)
//                if(fbuf->mst_info->masters[mst] == gl_my_rank)
//                    break;
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Changing data_sz from %d to %d", (int)data_sz, (int)(offt_range*(mst+1) - offset));
//            data_sz = offt_range*(mst+1) - offset;
//        }
//
//        new_chunk = new_write_chunk_rec(rank, offset, data_sz);
//
//        FARB_DBG(VERBOSE_ALL_LEVEL, "new chunk (r %d, off %llu, sz %llu)", new_chunk->rank, new_chunk->offset, new_chunk->data_sz);
//
//        if(dbitem->chunks == NULL){
//            dbitem->chunks = new_chunk;
//            dbitem->last = new_chunk;
//        } else {
//            write_chunk_rec_t *chunk;
//
//            chunk = dbitem->chunks;
//            while( (chunk != NULL) && (chunk->offset < new_chunk->offset))
//                chunk = chunk->next;
//
//            if(chunk == NULL){ //insert to the end
//                /*check for overlapping*/
//                if( (new_chunk->offset >= dbitem->last->offset) &&
//                    (new_chunk->offset < dbitem->last->offset + dbitem->last->data_sz)){
//
////                    FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: multiple writes of the same area (by processes %d and %d). \
////                            Overwrite order is not defined. ", dbitem->last->rank, new_chunk->rank);
//                    FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: write overlap (r%d, r%d): dbitem (%llu, %llu), new (%llu, %llu)",
//                             dbitem->last->rank, new_chunk->rank, dbitem->last->offset, dbitem->last->data_sz, new_chunk->offset, new_chunk->data_sz);
//                }
//                dbitem->last->next = new_chunk;
//                new_chunk->prev = dbitem->last;
//                dbitem->last = new_chunk;
//            } else {
//                    /*check for overlapping*/
//                    if(chunk->prev != NULL){
//                        if((new_chunk->offset >= chunk->prev->offset) &&
//                           (new_chunk->offset < chunk->prev->offset+chunk->prev->data_sz)){
////                            FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: multiple writes of the same area (by processes %d and %d). \
////                                    Overwrite order is not defined. ", chunk->prev->rank, new_chunk->rank);
//                            FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: write overlap (ranks %d,%d): dbitem (%llu, %llu), new (%llu, %llu)",
//                             chunk->prev->rank, new_chunk->rank, chunk->prev->offset, chunk->prev->data_sz, new_chunk->offset, new_chunk->data_sz);
//                        }
//                    }
//
//                   if( (new_chunk->offset + new_chunk->data_sz > chunk->offset) &&
//                       (new_chunk->offset + new_chunk->data_sz <= chunk->offset+chunk->data_sz)){
////                        FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: multiple writes of the same area (by processes %d and %d). \
////                                    Overwrite order is not defined. ", chunk->rank, new_chunk->rank);
//                        FARB_DBG(VERBOSE_ALL_LEVEL, "FARB Warning: write overlap (ranks %d, r%d): dbitem (%llu, %llu), new (%llu, %llu)",
//                             chunk->rank, new_chunk->rank, chunk->offset, chunk->data_sz, new_chunk->offset, new_chunk->data_sz);
//                    }
//
//                   //insert before this chunk
//                   if(chunk == dbitem->chunks){
//                        new_chunk->next = chunk;
//                        chunk->prev = new_chunk;
//                        dbitem->chunks = new_chunk;
//                   } else {
//                        chunk->prev->next = new_chunk;
//                        new_chunk->prev = chunk->prev;
//                        new_chunk->next = chunk;
//                        chunk->prev = new_chunk;
//                   }
//            }
//        }
//        fbuf->mst_info->iodb->updated_flag = 1;
//
//        }
//
//    FARB_DBG(VERBOSE_DBG_LEVEL, "List after update");
//    print_write_dbitem(dbitem);
//}
/*
    update the database with io request from reader
*/
//static void parse_read_ioreq(file_buffer_t *fbuf, farb_var_t *var, void *buf, int bufsz, int rank, MPI_Comm comm)
//{
//    read_db_item_t *dbitem;
//    read_chunk_rec_t *chunk, *tmp;
//
//    int mst1, mst2, mst;
//    MPI_Offset offset, data_sz, offt_range;
//    unsigned char *chbuf = (unsigned char*)buf;
//    size_t offt = 0;
//
//    /*Find corresponding record in the database*/
//    dbitem = fbuf->mst_info->iodb->ritems;
//    while(dbitem != NULL){
//        if( (dbitem->rank == rank) && (dbitem->comm == comm) )
//            break;
//        dbitem = dbitem->next;
//    }
//
//    if(dbitem == NULL){
//        dbitem = (read_db_item_t*)farb_malloc(sizeof(read_db_item_t));
//        assert(dbitem != NULL);
//        dbitem->rank = rank;
//        dbitem->comm = comm;
//        dbitem->next = NULL;
//        dbitem->prev = NULL;
//        dbitem->chunks = NULL;
//        dbitem->nchunks = 0;
//
//        //enqueue
//        if(fbuf->mst_info->iodb->ritems == NULL)
//            fbuf->mst_info->iodb->ritems = dbitem;
//        else{
//            dbitem->next = fbuf->mst_info->iodb->ritems;
//            fbuf->mst_info->iodb->ritems->prev = dbitem;
//            fbuf->mst_info->iodb->ritems = dbitem;
//        }
//        fbuf->mst_info->iodb->nritems++;
//    }
//    FARB_DBG(VERBOSE_DBG_LEVEL, "ritems %d",(int)fbuf->mst_info->iodb->nritems);
//
//    /*add mem chunks to the list for future matching.
//      Pick only those mem chunks whose offsets are within
//      the range I'm responsible for*/
//
//    if(has_unlim_dim(var)){
//        offt_range = 0;
//    } else {
//        int el_sz;
//        MPI_Type_size(var->dtype, &el_sz);
//        offt_range = (MPI_Offset)((last_1d_index(var->ndims, var->shape)+1)/fbuf->mst_info->nmasters);
//        offt_range *= el_sz;
//    }
//
//    FARB_DBG(VERBOSE_ALL_LEVEL, "offt_range %d", (int)offt_range);
//
//    while(offt != (size_t)bufsz){
//
//        offset = *(MPI_Offset*)(chbuf+offt);
//        offt+=sizeof(MPI_Offset);
//        data_sz = *(MPI_Offset*)(chbuf+offt);
//        offt+=sizeof(MPI_Offset);
//        if(offt_range == 0){
//            mst1 = mst2 = 0;
//        } else {
//            mst1 = (int)(offset/offt_range);
//            if(mst1 >= fbuf->mst_info->nmasters)
//                mst1 = fbuf->mst_info->nmasters - 1;
//            mst2 = (int)( (offset+data_sz-1)/offt_range);
//            if(mst2 >= fbuf->mst_info->nmasters)
//                mst2 = fbuf->mst_info->nmasters - 1;
//        }
//        FARB_DBG(VERBOSE_ALL_LEVEL, "--> (%d, %d): mst1 %d mst2 %d", (int)offset, (int)data_sz, mst1, mst2);
//        if( (fbuf->mst_info->masters[mst1] != gl_my_rank) && (fbuf->mst_info->masters[mst2] != gl_my_rank))
//            /*This master is not responsible for this data*/
//            continue;
//        else if(fbuf->mst_info->masters[mst1] != gl_my_rank){
//            for(mst = 0; mst < fbuf->mst_info->nmasters; mst++)
//                if(fbuf->mst_info->masters[mst] == gl_my_rank)
//                    break;
//            assert(mst != fbuf->mst_info->nmasters);
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Changing offset from %llu to %llu, datasz from %llu to %llu",
//                     offset, offt_range*mst, data_sz, data_sz - (offt_range*mst - offset));
//            data_sz -= offt_range*mst - offset;
//            offset = offt_range*mst;
//        } else if(fbuf->mst_info->masters[mst2] != gl_my_rank){
//            for(mst = 0; mst < fbuf->mst_info->nmasters; mst++)
//                if(fbuf->mst_info->masters[mst] == gl_my_rank)
//                    break;
//            assert(mst != fbuf->mst_info->nmasters);
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Changing data_sz from %d to %d", (int)data_sz, (int)(offt_range*(mst+1) - offset));
//            data_sz = offt_range*(mst+1) - offset;
//        }
//
//        chunk = (read_chunk_rec_t*)farb_malloc(sizeof(read_chunk_rec_t));
//        assert(chunk != NULL);
//        chunk->var_id = var->id;
//        chunk->offset = offset;
//        chunk->data_sz = data_sz;
//        chunk->next = NULL;
//        chunk->prev = NULL;
//
//        if(dbitem->chunks == NULL)
//            dbitem->chunks = chunk;
//        else{
//            tmp = dbitem->chunks;
//            while(tmp->next != NULL)
//                tmp = tmp->next;
//            tmp->next = chunk;
//            chunk->prev = tmp;
//        }
//        dbitem->nchunks++;
//        fbuf->mst_info->iodb->updated_flag = 1;
//    }
//
//    assert(dbitem->nchunks > 0);
//    FARB_DBG(VERBOSE_DBG_LEVEL, "rreq added");
//    print_read_dbitem(dbitem);
//}

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

    io_req_t *ioreq = (io_req_t*)farb_malloc(sizeof(io_req_t));
    assert(ioreq != NULL);
    int el_sz;
    MPI_Type_size(dtype, &el_sz);

    if(ndims > 0){
        ioreq->user_buf_sz = el_sz * (last_1d_index(ndims, count) + 1);
        ioreq->start = (MPI_Offset*)farb_malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->start != NULL);
        memcpy((void*)ioreq->start, (void*)start, sizeof(MPI_Offset)*ndims);
        ioreq->count = (MPI_Offset*)farb_malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->count != NULL);
        memcpy((void*)ioreq->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else{
        ioreq->start = NULL;
        ioreq->count = NULL;
        ioreq->user_buf_sz = el_sz;
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, user bufsz %d", id, var_id, (int)ioreq->user_buf_sz);
//    if( (rw_flag == FARB_WRITE) && (buffered)){
    //TODO buffered ioreq
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Dup user buffer");
////        ioreq->user_buf = farb_malloc((size_t)ioreq->user_buf_sz);
////        assert(ioreq->user_buf != NULL);
////        memcpy(ioreq->user_buf, buf, (size_t)ioreq->user_buf_sz);
//    } else
        ioreq->user_buf = buf;
    ioreq->next = NULL;
    ioreq->sent_flag = 0;
    ioreq->mem_chunks = NULL; //will init later
    ioreq->id = id;
    ioreq->var_id = var_id;
    ioreq->nchunks = 0;
    ioreq->get_sz = 0;
    ioreq->prev = NULL;
    ioreq->dtype = dtype;
    ioreq->rw_flag = rw_flag;
    return ioreq;
}

void delete_ioreq(file_buffer_t *fbuf, io_req_t **ioreq)
{
    contig_mem_chunk_t *chunk;
    farb_var_t *var = find_var(fbuf->vars, (*ioreq)->var_id);
    assert(var != NULL);

    if(*ioreq == fbuf->ioreqs)
        fbuf->ioreqs = fbuf->ioreqs->next;
    if((*ioreq)->next != NULL)
        (*ioreq)->next->prev = (*ioreq)->prev;
    if( (*ioreq)->prev != NULL)
        (*ioreq)->prev->next = (*ioreq)->next;

    if( (*ioreq)->rw_flag == FARB_READ)
        fbuf->rreq_cnt--;
    else
        fbuf->wreq_cnt--;

    if( (*ioreq)->count != NULL )
        farb_free((*ioreq)->count, var->ndims*sizeof(MPI_Offset));
    if((*ioreq)->start != NULL)
        farb_free((*ioreq)->start, var->ndims*sizeof(MPI_Offset));

    chunk = (*ioreq)->mem_chunks;
    while(chunk != NULL){
        (*ioreq)->mem_chunks = chunk->next;
        farb_free(chunk, sizeof(contig_mem_chunk_t));
        chunk = (*ioreq)->mem_chunks;
    }
    farb_free((*ioreq), sizeof(io_req_t));
}


void send_ioreqs(file_buffer_t *fbuf)
{
    io_req_t *ioreq;
    contig_mem_chunk_t *chunk;
    unsigned char **sbuf;
    size_t *offt, *bufsz, *nchuncks_offt;
    int errno, i;
    MPI_Request *sreqs;
    int nmasters = fbuf->mst_info->nmasters;
    size_t dflt_sz = 256*sizeof(MPI_Offset);
    farb_var_t *var;
    MPI_Offset cur_offt, cur_dsz, tmp_dsz;
    int el_sz;
    size_t offt_range;
    int mst;
    int *mst_flag;
    unsigned int *nchunks;
    double t_start = MPI_Wtime();

    assert(nmasters > 0);
    /*init*/
    sbuf = (unsigned char**)farb_malloc(nmasters*sizeof(unsigned char*));
    assert(sbuf != NULL);
    sreqs = (MPI_Request*)farb_malloc(nmasters*sizeof(MPI_Request));
    assert(sreqs != NULL);
    offt = (size_t*)farb_malloc(nmasters*sizeof(size_t));
    assert(offt != NULL);
    bufsz = (size_t*)farb_malloc(nmasters*sizeof(size_t));
    assert(bufsz != NULL);
    nchuncks_offt = (size_t*)farb_malloc(nmasters*sizeof(size_t));
    assert(nchuncks_offt != NULL);
    mst_flag = (int*)farb_malloc(nmasters*sizeof(int));
    assert(mst_flag != NULL);
    nchunks = (unsigned int*)farb_malloc(nmasters*sizeof(unsigned int));
    assert(nchunks != NULL);

    for(i = 0; i < nmasters; i++){
        sreqs[i] = MPI_REQUEST_NULL;
        //offt[i] = 0;
        sbuf[i] = (unsigned char*)farb_malloc(dflt_sz);
        assert(sbuf[i] != NULL);
        *(MPI_Offset*)sbuf[i] = (MPI_Offset)fbuf->ncid;
        offt[i] = sizeof(MPI_Offset);
        bufsz[i] = dflt_sz;
        mst_flag[i] = 0;
    }

    /*Pack requests*/
    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            ioreq = ioreq->next;
            continue;
        }
        if(ioreq->rw_flag == FARB_READ)
            FARB_DBG(VERBOSE_DBG_LEVEL, "Pack read request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);
        else
            FARB_DBG(VERBOSE_DBG_LEVEL, "Pack write request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);

        /*Store varid, rw flag and remember the offset where number of chunks will
        be written
        */
        for(i = 0; i < nmasters; i++){
            /*extend the buf if needed*/
            if(offt[i]+sizeof(MPI_Offset)*3 > bufsz[i]){
                unsigned char* tmp = realloc(sbuf[i], bufsz[i]+dflt_sz);
                assert(tmp != NULL);
                sbuf[i] = tmp;
                FARB_DBG(VERBOSE_DBG_LEVEL, "buf reallocd fr %lu to %lu", bufsz[i], bufsz[i]+dflt_sz);
                bufsz[i] += dflt_sz;
                gl_conf.malloc_size += dflt_sz;

            }
            //FARB_DBG(VERBOSE_DBG_LEVEL, "offt %lu", offt[i]);
            //assert(sbuf[i] != NULL);
            *(MPI_Offset*)(sbuf[i] + offt[i]) = (MPI_Offset)ioreq->rw_flag;
            offt[i] += sizeof(MPI_Offset);
            *(MPI_Offset*)(sbuf[i] + offt[i]) = (MPI_Offset)ioreq->var_id;
            offt[i] += sizeof(MPI_Offset);
            /*remember where to save the written nchunks later*/
            nchuncks_offt[i] = offt[i];
            offt[i] += sizeof(MPI_Offset);
            nchunks[i] = 0; //reset
        }

        var = find_var(fbuf->vars, ioreq->var_id);
        /*Find out the offset range size for this var*/
        MPI_Type_size(var->dtype, &el_sz);
        assert(el_sz > 0);
        if(has_unlim_dim(var)){
            int nelems = 1;
            assert(var->shape[0] == FARB_UNLIMITED);
            for(i = 1; i < var->ndims; i++)
                nelems *= var->shape[i];
            offt_range = nelems * el_sz * UNLIM_NELEMS_RANGE;
        } else {
            //FARB_DBG(VERBOSE_ALL_LEVEL, "last idx %d, offt %d", (int)last_1d_index(var->ndims, var->shape), (int)(last_1d_index(var->ndims, var->shape)*el_sz));
            offt_range = (MPI_Offset)((last_1d_index(var->ndims, var->shape) + 1)/fbuf->mst_info->nmasters);
            offt_range *= el_sz;
        }
        FARB_DBG(VERBOSE_DBG_LEVEL, "offt_range %d", (int)offt_range);
        assert(offt_range > 0);

        chunk = ioreq->mem_chunks;
        while(chunk != NULL){

            cur_offt = chunk->offset;
            cur_dsz = chunk->data_sz;

            int tmp;
            while(cur_dsz){
                /*found out the master who's responsible for
                 this data.*/
                tmp = (int)(cur_offt/offt_range);
                mst = tmp%nmasters;
                /*if the last written bit for this chunk
                  falls inside the range of the next master
                  we have to split the current chunk*/
                if(cur_offt+cur_dsz-1 >= (tmp+1)*offt_range){
                    tmp_dsz = (tmp+1)*offt_range - cur_offt;
                } else
                    tmp_dsz = cur_dsz;

                /*extend the buf if needed*/
                if(offt[mst]+sizeof(MPI_Offset)*2 > bufsz[mst]){
                    unsigned char* tmp = realloc(sbuf[mst], bufsz[mst]+dflt_sz);
                    assert(tmp != NULL);
                    sbuf[mst] = tmp;
                    FARB_DBG(VERBOSE_DBG_LEVEL, "buf reallocd fr %lu to %lu", bufsz[mst], bufsz[mst]+dflt_sz);
                    bufsz[mst] += dflt_sz;
                    gl_conf.malloc_size += dflt_sz;
                }
                FARB_DBG(VERBOSE_ALL_LEVEL, "  chunk (r %d, off %lld, sz %lld)", fbuf->mst_info->masters[mst],
                         cur_offt, tmp_dsz);
                *(MPI_Offset*)(sbuf[mst]+offt[mst]) = cur_offt;
                offt[mst] += sizeof(MPI_Offset);
                *(MPI_Offset*)(sbuf[mst]+offt[mst]) = tmp_dsz;
                offt[mst] += sizeof(MPI_Offset);
                nchunks[mst]++;
                mst_flag[mst] = 1;
                cur_offt += tmp_dsz;
                cur_dsz -= tmp_dsz;
            }

            chunk = chunk->next;
        }
        /*write number of chunks stored*/
        for(i = 0; i < nmasters; i++){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Stored %u chunks for mst %d", nchunks[i], fbuf->mst_info->masters[i]);
            *(MPI_Offset*)(sbuf[i] + nchuncks_offt[i]) = (MPI_Offset)nchunks[i];
        }
        ioreq->sent_flag = 1;
        ioreq = ioreq->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Before sending reqs");
    /*Send the request to the master*/
    for(i = 0; i < nmasters; i++){
         if(!mst_flag[i])
            continue;   //nothing to send
        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[i] == gl_my_rank))
            parse_ioreqs(sbuf[i], (int)offt[i], gl_my_rank, gl_comps[gl_my_comp_id].comm);
        else{
            assert(gl_comps[fbuf->writer_id].comm != MPI_COMM_NULL);
            FARB_DBG(VERBOSE_DBG_LEVEL, "Send reqs to mst %d (bufsz %lu)", fbuf->mst_info->masters[i], offt[i]);
            errno = MPI_Isend((void*)sbuf[i], (int)offt[i], MPI_BYTE, fbuf->mst_info->masters[i], IO_REQS_TAG, gl_comps[fbuf->writer_id].comm, &sreqs[i]);
            CHECK_MPI(errno);
            gl_stats.nmatching_msg_sent++;
            gl_stats.accum_msg_sz += bufsz[i];
        }
    }
    errno = MPI_Waitall(nmasters, sreqs, MPI_STATUSES_IGNORE);
    CHECK_MPI(errno);

    FARB_DBG(VERBOSE_DBG_LEVEL, "Sent reqs");


    /*dealloc*/
    farb_free(sreqs, fbuf->mst_info->nmasters*sizeof(MPI_Request));
    for(i = 0; i < nmasters; i++){
        farb_free(sbuf[i], bufsz[i]);
    }
    farb_free(sbuf, nmasters*sizeof(unsigned char*));
    farb_free(offt, nmasters*sizeof(size_t));
    farb_free(bufsz, nmasters*sizeof(size_t));
    farb_free(nchuncks_offt, nmasters*sizeof(size_t));
    farb_free(mst_flag, sizeof(int) * nmasters);
    farb_free(nchunks, sizeof(unsigned int)*nmasters);

    FARB_DBG(VERBOSE_DBG_LEVEL, "Stat: Time to pack and send reqs %.4f", MPI_Wtime() - t_start);
}

//void send_ioreqs(file_buffer_t *fbuf)
//{
//    io_req_t *ioreq;
//    contig_mem_chunk_t *chunk;
//    unsigned char *sbuf = NULL;
//    size_t sz = 0, offt = 0;
//    int errno, i;
//   // int *mst_flag, mst1, mst2, i, errno;
//    MPI_Request *sreqs;
// //   MPI_Offset offt_range;
////    farb_var_t *var;
//
////    mst_flag = (int*)farb_malloc(fbuf->mst_info->nmasters*sizeof(int));
////    assert(mst_flag != NULL);
//    sreqs = (MPI_Request*)farb_malloc(fbuf->mst_info->nmasters*sizeof(MPI_Request));
//    assert(sreqs != NULL);
//    for(i = 0; i < fbuf->mst_info->nmasters; i++){
// //       mst_flag[i] = 0;
//        sreqs[i] = MPI_REQUEST_NULL;
//    }
//
//    /*first find out the size of the buffer*/
//    sz = sizeof(MPI_Offset); //ncid
//    ioreq = fbuf->ioreqs;
//    while(ioreq != NULL){
//        if(ioreq->sent_flag){
//            ioreq = ioreq->next;
//            continue;
//        }
//        sz += sizeof(MPI_Offset)*3 + sizeof(MPI_Offset)*2*ioreq->nchunks; //var_id+rw_flag+nchunks+chunks
//        ioreq = ioreq->next;
//    }
//
//    sbuf = farb_malloc(sz);
//    assert(sbuf != NULL);
//
//    FARB_DBG(VERBOSE_DBG_LEVEL, "Bufsz %d", (int)sz);
//
//    offt = 0;
//
//    *(MPI_Offset*)sbuf = (MPI_Offset)fbuf->ncid;
//    offt += sizeof(MPI_Offset);
//
//    /*Pack requests*/
//    ioreq = fbuf->ioreqs;
//    while(ioreq != NULL){
//        if(ioreq->sent_flag){
//            ioreq = ioreq->next;
//            continue;
//        }
//        if(ioreq->rw_flag == FARB_READ)
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Pack read request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);
//        else
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Pack write request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);
//
//        *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->rw_flag;
//        offt += sizeof(MPI_Offset);
//        *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->var_id;
//        offt += sizeof(MPI_Offset);
//        *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->nchunks;
//        offt += sizeof(MPI_Offset);
//
//        chunk = ioreq->mem_chunks;
//        while(chunk != NULL){
//            *(MPI_Offset*)(sbuf+offt) = chunk->offset;
//            offt += sizeof(MPI_Offset);
//            *(MPI_Offset*)(sbuf+offt) = chunk->data_sz;
//            offt += sizeof(MPI_Offset);
//            chunk = chunk->next;
//        }
//        ioreq->sent_flag = 1;
//        ioreq = ioreq->next;
//    }
//
//    assert(offt == sz);
//    FARB_DBG(VERBOSE_DBG_LEVEL, "Before sending reqs");
//    /*Send the request to the master*/
//    for(i = 0; i < fbuf->mst_info->nmasters; i++){
//         errno = MPI_SUCCESS;
//
//        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[i] == gl_my_rank))
//            parse_ioreqs(sbuf, (int)sz, gl_my_rank, gl_comps[gl_my_comp_id].comm);
//        else{
//            assert(gl_comps[fbuf->writer_id].comm != MPI_COMM_NULL);
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Send reqs to mst %d", fbuf->mst_info->masters[i]);
//            errno = MPI_Isend((void*)sbuf, (int)sz, MPI_BYTE, fbuf->mst_info->masters[i], IO_REQS_TAG, gl_comps[fbuf->writer_id].comm, &sreqs[i]);
//            CHECK_MPI(errno);
//            gl_stats.nmatching_msg_sent++;
//            gl_stats.accum_msg_sz += sz;
//        }
//    }
//    errno = MPI_Waitall(fbuf->mst_info->nmasters, sreqs, MPI_STATUSES_IGNORE);
//    CHECK_MPI(errno);
//
//    FARB_DBG(VERBOSE_DBG_LEVEL, "Sent reqs");
////    while(!completed){
////        errno = MPI_Testall(fbuf->mst_info->nmasters, sreqs, &completed, MPI_STATUSES_IGNORE);
////        CHECK_MPI(errno);
////        progress_io_matching();
////    }
//
//    farb_free(sbuf, sz);
//    //farb_free(mst_flag, fbuf->mst_info->nmasters*sizeof(int));
//    farb_free(sreqs, fbuf->mst_info->nmasters*sizeof(MPI_Request));
//}


/* Reader ranks and writer ranks send their read/write
   io request to master(s)
*/
//void send_ioreq(file_buffer_t *fbuf, io_req_t *ioreq)
//{
//    contig_mem_chunk_t *chunk;
//    unsigned char *sbuf = NULL;
//    size_t sz = 0, offt = 0;
//    int *mst_flag, mst1, mst2, i, errno;
//    MPI_Request *sreqs;
//    MPI_Offset offt_range;
//    farb_var_t *var;
//    if(ioreq->sent_flag)
//        return;
//    if(ioreq->rw_flag == FARB_READ)
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Send read request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);
//    else
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Send write request to master for var %d (%s)", ioreq->var_id, fbuf->file_path);
//
//    mst_flag = (int*)farb_malloc(fbuf->mst_info->nmasters*sizeof(int));
//    assert(mst_flag != NULL);
//    sreqs = (MPI_Request*)farb_malloc(fbuf->mst_info->nmasters*sizeof(MPI_Request));
//    assert(sreqs != NULL);
//    for(i = 0; i < fbuf->mst_info->nmasters; i++){
//        mst_flag[i] = 0;
//        sreqs[i] = MPI_REQUEST_NULL;
//    }
//
//    /*Pack the read request.*/
//    offt = 0;
//    //sz = sizeof(ncid)+sizeof(ioreq->var_id)+sizeof(MPI_Offset)*2*ioreq->nchunks;
//    sz = sizeof(MPI_Offset)*3+sizeof(MPI_Offset)*2*ioreq->nchunks;
//    sbuf = farb_malloc(sz);
//    assert(sbuf != NULL);
//
//    FARB_DBG(VERBOSE_ALL_LEVEL, "Bufsz %d", (int)sz);
//
//    *(MPI_Offset*)sbuf = (MPI_Offset)fbuf->ncid;
//    offt += sizeof(MPI_Offset);
//    *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->rw_flag;
//    offt += sizeof(MPI_Offset);
//    *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->var_id;
//    offt += sizeof(MPI_Offset);
//
//    var = find_var(fbuf->vars, ioreq->var_id);
//    assert(var != NULL);
//
//    for(i = 0; i < var->ndims; i++)
//        FARB_DBG(VERBOSE_DBG_LEVEL, "---> start %d, count %d", (int)ioreq->start[i], (int)ioreq->count[i]);
//
//    if(has_unlim_dim(var)){
//        /*For vars with unlimited dimension all the reqs
//          will be sent to master 0*/
//        offt_range = 0;
//        FARB_DBG(VERBOSE_DBG_LEVEL, "Req for a var with unlim dim");
//    } else {
//        int el_sz;
//        MPI_Type_size(var->dtype, &el_sz);
//        offt_range = (MPI_Offset)(( last_1d_index(var->ndims, var->shape) + 1 )/fbuf->mst_info->nmasters);
//        offt_range *= el_sz;
//    }
//    FARB_DBG(VERBOSE_ALL_LEVEL, "Offset range is %d", (int)offt_range);
//
//    chunk = ioreq->mem_chunks;
//    while(chunk != NULL){
//        *(MPI_Offset*)(sbuf+offt) = chunk->offset;
//        offt += sizeof(MPI_Offset);
//        *(MPI_Offset*)(sbuf+offt) = chunk->data_sz;
//        offt += sizeof(MPI_Offset);
//
//        /*Set flag for each master to whom will need to send this request*/
//        if(offt_range == 0)
//            mst1 = mst2 = 0;
//        else {
//            mst1 = (int)(chunk->offset/offt_range);
//            if(mst1 >= fbuf->mst_info->nmasters)
//                mst1 = fbuf->mst_info->nmasters - 1;
//            mst2 = (int)((chunk->offset + chunk->data_sz - 1)/offt_range);
//            if(mst2 >= fbuf->mst_info->nmasters)
//                mst2 = fbuf->mst_info->nmasters - 1;
//        }
//        FARB_DBG(VERBOSE_ALL_LEVEL, "-> (%d, %d): mst1 %d mst2 %d", (int)chunk->offset, (int)chunk->data_sz, mst1, mst2);
//        for(i = mst1; i <= mst2; i++)
//            mst_flag[i] = 1;
//        chunk = chunk->next;
//    }
//    assert(offt == sz);
//
//    /*Send the request to the master*/
//    for(i = 0; i < fbuf->mst_info->nmasters; i++){
//        if(!mst_flag[i])
//            continue;
//        errno = MPI_SUCCESS;
//
//        if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->mst_info->masters[i] == gl_my_rank))
//            parse_ioreq(sbuf, (int)sz, gl_my_rank, gl_comps[gl_my_comp_id].comm);
//        else{
//            assert(gl_comps[fbuf->writer_id].comm != MPI_COMM_NULL);
//            FARB_DBG(VERBOSE_DBG_LEVEL, "Send req %u to mst %d", ioreq->id, fbuf->mst_info->masters[i]);
//            errno = MPI_Isend((void*)sbuf, (int)sz, MPI_BYTE, fbuf->mst_info->masters[i], IO_REQ_TAG, gl_comps[fbuf->writer_id].comm, &sreqs[i]);
//            CHECK_MPI(errno);
//            gl_stats.nmatching_msg_sent++;
//            gl_stats.accum_msg_sz += sz;
//        }
//    }
//    errno = MPI_Waitall(fbuf->mst_info->nmasters, sreqs, MPI_STATUSES_IGNORE);
//    CHECK_MPI(errno);
//    ioreq->sent_flag = 1;
////    while(!completed){
////        errno = MPI_Testall(fbuf->mst_info->nmasters, sreqs, &completed, MPI_STATUSES_IGNORE);
////        CHECK_MPI(errno);
////        progress_io_matching();
////    }
//
//    farb_free(sbuf, sz);
//    farb_free(mst_flag, fbuf->mst_info->nmasters*sizeof(int));
//    farb_free(sreqs, fbuf->mst_info->nmasters*sizeof(MPI_Request));
//}

void match_ioreqs_all(int rw_flag)
{
    file_buffer_t *fbuf;
    io_req_t *ioreq;
    int fbuf_cnt = 0, tmp_cnt;

    /*Note: rw_flag should be FARB_WRITE (checked this in higher level)*/
    assert(rw_flag == FARB_WRITE);

    FARB_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs all");

    /*Do preparations*/
    fbuf = gl_filebuf_list;
    while(fbuf != NULL){
        if( (fbuf->iomode != FARB_IO_MODE_MEMORY) ||
            (!fbuf->explicit_match)){
            // || ( (rw_flag == FARB_WRITE) && (fbuf->writer_id != gl_my_comp_id))){

            fbuf = fbuf->next;
            continue;
        }
        /*Count for how many files we need to match ioreqs.*/
        fbuf_cnt++;
        /*Init master's db*/
        if(fbuf->mst_info->is_master_flag){
            /*Check there is no unfinished matching process*/
            if(fbuf->is_matching_flag == 1){
                FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Error: Matching for %s has not completed yet. Cannot do match all.", fbuf->file_path);
                assert(fbuf->is_matching_flag == 0);
            }
        }
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will match for %s", fbuf->file_path);

        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if(ioreq->rw_flag==FARB_READ && fbuf->writer_id==gl_my_comp_id){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: cannot call farb_match_ioreqs_all when there are uncompleted intra-component I/O requests ");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
            ioreq = ioreq->next;
        }
        /*Send ioreqs to master(s)*/
        send_ioreqs(fbuf);
        fbuf->done_matching_flag = 0;
        fbuf->is_matching_flag = 1;
        fbuf = fbuf->next;
    }

    if(fbuf_cnt == 0){
        FARB_DBG(VERBOSE_DBG_LEVEL, "There are no files to match ioreqs.");
        return;
    } else
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will match ioreqs for %d files", fbuf_cnt);
    int counter = 0;
    /*Keep iterating over files and trying to progress I/O*/
    while(fbuf_cnt){
        tmp_cnt = 0;
        fbuf = gl_filebuf_list;
        while(fbuf != NULL){
            if( (fbuf->iomode != FARB_IO_MODE_MEMORY) ||
                (!fbuf->explicit_match)){
                //|| ( (rw_flag == FARB_WRITE) && (fbuf->writer_id != gl_my_comp_id))){

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
                counter++;
                if(counter % 300 == 0){
                    do_matching(fbuf, 0);
                    counter = 0;
                }
            }

            fbuf = fbuf->next;
        }
        if(tmp_cnt == fbuf_cnt)
            break;
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished matching I/O for all files");

    //reset all flags
    fbuf = gl_filebuf_list;
    while(fbuf != NULL){
        if( (fbuf->iomode != FARB_IO_MODE_MEMORY) ||(!fbuf->explicit_match) ){
            //|| ( (rw_flag == FARB_WRITE) && (fbuf->writer_id != gl_my_comp_id))){
            fbuf = fbuf->next;
            continue;
        }

        fbuf->is_matching_flag = 0;
        fbuf->mst_info->nrranks_completed = 0;
        assert(fbuf->mst_info->nwranks_completed == 0);
        fbuf->done_matching_flag = 0;
        fbuf = fbuf->next;
    }
}

int match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag)
{
    double t_start;

    FARB_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs for file %d, intracomp %d", fbuf->ncid, intracomp_io_flag);

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
        FARB_DBG(VERBOSE_DBG_LEVEL, "PFARB Warning: ps has no requests for file %d", fbuf->ncid);
        if(fbuf->reader_id == gl_my_comp_id){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Completed all rreqs. Notify master.");
            /*Notify master 0 that all my read io
            requests for this file have been completed*/
            int errno = MPI_Send(&(fbuf->ncid), 1, MPI_INT, fbuf->root_writer, READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
            CHECK_MPI(errno);
            gl_stats.nmsg_sent++;
            return 0;
        }
    } else{
        FARB_DBG(VERBOSE_DBG_LEVEL, "Total %d rreqs and %d wreqs to match", fbuf->rreq_cnt, fbuf->wreq_cnt);
        send_ioreqs(fbuf);
//        ioreq = fbuf->ioreqs;
//        while(ioreq != NULL){
//            /*Forward the info about the request to writer's master rank(s)*/
//            send_ioreq(fbuf, ioreq);
//            ioreq = ioreq->next;
//        }
    }
    //reset
    fbuf->done_matching_flag = 0;
    /*Have to check for both things otherwise writers sometime hang
    (maybe they receive IO_CLOSE_FILE_TAG before READ_DONE_TAG by mistake?)*/
    int counter = 0;
    while(!fbuf->done_matching_flag){
        progress_io_matching();


        if( (fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag  ){
            counter++;
            if(counter % 300 == 0){
                do_matching(fbuf, intracomp_io_flag);
                counter = 0;
            }
        }
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);
    //reset
    if( (fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag ){
        if(fbuf->mst_info->nrranks_completed != 0){
            assert(fbuf->mst_info->nrranks_completed == fbuf->mst_info->nrranks_opened);
            fbuf->mst_info->nrranks_completed = 0;
        }
        if(fbuf->mst_info->nwranks_completed != 0){
            FARB_DBG(VERBOSE_DBG_LEVEL, "compl %d, open %d",fbuf->mst_info->nwranks_completed, fbuf->mst_info->nwranks_opened );
            assert(fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened);
            fbuf->mst_info->nwranks_completed = 0;
        }
        fbuf->done_matching_flag = 0;
    }

    /*Make readers synch before they complete. Otherwise it may happen
      that if readers call match_ioreqs multiple times, the ioreqs
      from two different match processes may mix up on the writer's side.
      It may lead to a dead lock situation.*/

      //if(fbuf->reader_id == gl_my_comp_id){
      FARB_DBG(VERBOSE_DBG_LEVEL, "Barrier before completing matching for %s", fbuf->file_path);
      MPI_Barrier(fbuf->comm);
      //}
    gl_stats.accum_match_time += MPI_Wtime() - t_start;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Stat: Time for matching %.4f", MPI_Wtime() - t_start);
    return 0;
}

/*writer->reader*/
static void send_data_wrt2rdr(void* buf, int bufsz)
{
    int ncid, var_id, rdr_rank, errno, i;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    farb_var_t *var;
    size_t rofft = 0, sofft = 0;
    unsigned char *sbuf = NULL;
    size_t sbufsz = 0;
    MPI_Offset read_offt, read_sz, dbuf_offt;
    contig_mem_chunk_t *chunk;
    int def_el_sz, req_el_sz, intracomp_flag;
    unsigned char *rbuf = (unsigned char*)buf;

    rdr_rank = (int)(*(MPI_Offset*)rbuf);
    rofft += sizeof(MPI_Offset);
    intracomp_flag = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Sending data to rank %d (intracomp flag %d)", rdr_rank, intracomp_flag);
    ncid = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    assert(fbuf != NULL);

    /*First compute the size of buf to allocate*/
    //skip var_id and start_offt
    rofft += sizeof(MPI_Offset)*2;
    sbufsz = sizeof(MPI_Offset); //ncid
    while(rofft < bufsz){
        read_sz = *(MPI_Offset*)(rbuf+rofft);
        //padding
        read_sz += read_sz%sizeof(MPI_Offset);
        //var_id+start_offt+chunk_sz+chunk
        sbufsz += sizeof(MPI_Offset)+sizeof(MPI_Offset)*2+(size_t)read_sz;
        //skip until next chunk_sz
        rofft += sizeof(MPI_Offset)+sizeof(MPI_Offset)*2;
    }

    FARB_DBG(VERBOSE_ALL_LEVEL, "Alloc buf of sz %d", (int)sbufsz);
    sbuf = farb_malloc(sbufsz);
    assert(sbuf != NULL);
    *(MPI_Offset*)sbuf = (MPI_Offset)ncid;
    sofft = sizeof(MPI_Offset);

    //rewind
    rofft = sizeof(MPI_Offset)*3;
    while(rofft != bufsz){
        var_id = (int)(*(MPI_Offset*)(rbuf+rofft));
        rofft += sizeof(MPI_Offset);

        read_offt = *(MPI_Offset*)(rbuf+rofft);
        rofft += sizeof(MPI_Offset);
        read_sz = *(MPI_Offset*)(rbuf+rofft);
        rofft += sizeof(MPI_Offset);

        /*Save info about the memory chunk*/
        *(MPI_Offset*)(sbuf+sofft) = (MPI_Offset)var_id;
        sofft += sizeof(MPI_Offset);
        *(MPI_Offset*)(sbuf+sofft) = read_offt;
        sofft += sizeof(MPI_Offset);
        *(MPI_Offset*)(sbuf+sofft) = read_sz;
        sofft += sizeof(MPI_Offset);

        FARB_DBG(VERBOSE_DBG_LEVEL, "Will copy (var %d, offt %d, sz %d), padding %d", var_id, (int)read_offt, (int)read_sz, (int)(read_sz%sizeof(MPI_Offset)));

        /*Find the ioreq that has info about this mem chunk*/
        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if( (ioreq->var_id == var_id) && (ioreq->rw_flag == FARB_WRITE)){
                chunk = ioreq->mem_chunks;
                while(chunk != NULL){
                    if( (read_offt >= chunk->offset) &&
                        (read_offt < chunk->offset+chunk->data_sz) &&
                        (read_offt+read_sz <= chunk->offset+chunk->data_sz))
                        break;

                    chunk = chunk->next;
                }
                if(chunk != NULL)
                    break;
            }
            ioreq = ioreq->next;
        }
        //FARB_DBG(VERBOSE_ALL_LEVEL, "BEfore assert");
        assert(ioreq != NULL);
        //FARB_DBG(VERBOSE_ALL_LEVEL, "After assert");
        assert(chunk != NULL);

        var = find_var(fbuf->vars, var_id);
        /*Copy data: if the var datatype passed in I/O function differs
        from the datatype with which the var was defined then we will
        need to do type conversion*/
        MPI_Type_size(ioreq->dtype, &req_el_sz);
        MPI_Type_size(var->dtype, &def_el_sz);
        if(req_el_sz == def_el_sz){
            assert(read_offt - ioreq->mem_chunks->offset >= 0);
            dbuf_offt = chunk->usrbuf_offset + read_offt - chunk->offset;
            //FARB_DBG(VERBOSE_ALL_LEVEL, "dbuf offt %d, chunk usr buf offt %d (usrbuf sz %d)", (int)dbuf_offt, (int)chunk->usrbuf_offset, (int)ioreq->user_buf_sz);
            assert(dbuf_offt < ioreq->user_buf_sz);
            memcpy(sbuf+sofft, (unsigned char*)ioreq->user_buf+dbuf_offt, (size_t)read_sz);
            sofft += (size_t)read_sz+/*padding*/(size_t)(read_sz%sizeof(MPI_Offset));
        } else {
            void *cpbuf;
            MPI_Offset offt_within_chunk;
            int nelems_to_read;

            FARB_DBG(VERBOSE_ALL_LEVEL, "Converting from %d-bit to %d bit type for var %d", req_el_sz, def_el_sz, var_id);
            assert(read_offt - ioreq->mem_chunks->offset >= 0);
            assert((read_offt - chunk->offset)%def_el_sz == 0);
            offt_within_chunk = (MPI_Offset)((read_offt - chunk->offset)/def_el_sz)*req_el_sz;
            nelems_to_read = (int)(read_sz/def_el_sz);

            cpbuf = farb_malloc(read_sz);
            assert(cpbuf != NULL);
            dbuf_offt = chunk->usrbuf_offset + offt_within_chunk;

            if(var->dtype == MPI_FLOAT){
                double *ubuf;
                /*MPI_DOUBLE->MPI_FLOAT*/
                assert(ioreq->dtype == MPI_DOUBLE);
                ubuf = (double*)((unsigned char*)ioreq->user_buf+dbuf_offt);
                for(i = 0; i < nelems_to_read; i++){
                    ((float*)cpbuf)[i] = (float)(ubuf[i]);
                    FARB_DBG(VERBOSE_ALL_LEVEL, "Send: uval %.10f, bval %.10f",ubuf[i], ((float*)cpbuf)[i]);
                }
            } else if (var->dtype == MPI_DOUBLE){
                float *ubuf;
              //  FARB_DBG(VERBOSE_DBG_LEVEL, "req_el_sz %d, def_el_sz %d, req dt %d, def dt %d", req_el_sz, def_el_sz, (int)ioreq->dtype, (int)var->dtype);
                /*MPI_FLOAT->MPI_DOUBLE*/
                assert(ioreq->dtype == MPI_FLOAT);
                ubuf = (float*)((unsigned char*)ioreq->user_buf+dbuf_offt);
                for(i = 0; i < nelems_to_read; i++){
                    ((double*)cpbuf)[i] = (double)(ubuf[i]);
                    FARB_DBG(VERBOSE_ALL_LEVEL, "Send: uval %.10f, bval %.10f",ubuf[i], ((double*)cpbuf)[i]);
                }
            } else {
                FARB_DBG(VERBOSE_ERROR_LEVEL, "This conversion type is not supported. Aborting.");
                assert(0);
            }
            memcpy(sbuf+sofft, cpbuf, (size_t)read_sz);
            sofft += (size_t)read_sz+/*padding*/(size_t)(read_sz%sizeof(MPI_Offset));
            farb_free(cpbuf, read_sz);
        }
    }
    FARB_DBG(VERBOSE_ALL_LEVEL, "Sofft %d, sbufsz %d", (int)sofft, (int)sbufsz);
    assert(sofft == sbufsz);
    if(intracomp_flag)
        errno = MPI_Send((void*)sbuf, (int)sbufsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[gl_my_comp_id].comm);
    else
        errno = MPI_Send((void*)sbuf, (int)sbufsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm);
    CHECK_MPI(errno);
    gl_stats.nmatching_msg_sent++;
    gl_stats.accum_msg_sz += sbufsz;
    farb_free(sbuf, sbufsz);
}

/*writer->reader*/
static void recv_data_rdr(void* buf, int bufsz)
{
    int ncid, var_id, i;
    int def_el_sz, req_el_sz;
    farb_var_t *var;
    io_req_t *ioreq = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    MPI_Offset read_offt, read_sz, dbuf_offt;
    contig_mem_chunk_t *chunk;
    unsigned char *chbuf = (unsigned char*)buf;
    ncid = (int)(*(MPI_Offset*)(chbuf));
    offt += sizeof(MPI_Offset);

    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    assert(fbuf != NULL);

    FARB_DBG(VERBOSE_DBG_LEVEL, "Received data for file %s (ncid %d)", fbuf->file_path, fbuf->ncid);
    while(offt != bufsz){
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        read_offt = *(MPI_Offset*)(chbuf+offt);
        offt += sizeof(MPI_Offset);
        read_sz = *(MPI_Offset*)(chbuf+offt);
        offt += sizeof(MPI_Offset);

        FARB_DBG(VERBOSE_DBG_LEVEL, "Will recv (var %d, offt %d, sz %d)", var_id, (int)read_offt, (int)read_sz);

        /*Find the ioreq*/
        ioreq = fbuf->ioreqs;
        while(ioreq != NULL){
            if(ioreq->var_id == var_id && ioreq->rw_flag == FARB_READ){
                chunk = ioreq->mem_chunks;
                while(chunk != NULL){
                    if( (read_offt >= chunk->offset) &&
                        (read_offt < chunk->offset+chunk->data_sz) &&
                        (read_offt+read_sz <= chunk->offset+chunk->data_sz))
                        break;

                    chunk = chunk->next;
                }
                if(chunk != NULL)
                    break;
            }
            ioreq = ioreq->next;
        }
        assert(ioreq != NULL);
        assert(chunk != NULL);

        var = find_var(fbuf->vars, var_id);
        /*Copy data*/
        MPI_Type_size(ioreq->dtype, &req_el_sz);
        MPI_Type_size(var->dtype, &def_el_sz);
        if(req_el_sz == def_el_sz){
            assert(read_offt - ioreq->mem_chunks->offset >= 0);
            dbuf_offt = chunk->usrbuf_offset + read_offt - chunk->offset;
            FARB_DBG(VERBOSE_ALL_LEVEL, "dbuf offt %d", (int)dbuf_offt);
            assert(dbuf_offt < ioreq->user_buf_sz);
            memcpy((unsigned char*)(ioreq->user_buf)+dbuf_offt, (unsigned char*)buf+offt, (size_t)read_sz);
            offt += (size_t)read_sz+/*padding*/(size_t)(read_sz%sizeof(MPI_Offset));
            ioreq->get_sz += read_sz;
        } else {
            void *cpbuf;
            MPI_Offset offt_within_chunk;
            int nelems_to_read;

            FARB_DBG(VERBOSE_ALL_LEVEL, "Converting from %d-bit to %d bit type for var %d", req_el_sz, def_el_sz, var_id);
            assert(read_offt - ioreq->mem_chunks->offset >= 0);
            assert((read_offt - chunk->offset)%def_el_sz == 0);
            offt_within_chunk = (MPI_Offset)((read_offt - chunk->offset)/def_el_sz)*req_el_sz;
            dbuf_offt = chunk->usrbuf_offset + offt_within_chunk;

            cpbuf = farb_malloc(read_sz);
            assert(cpbuf != NULL);

            /*Copy the data to temporary buffer and then do the conversion*/
            memcpy(cpbuf, (unsigned char*)buf+offt, (size_t)read_sz);
            offt += (size_t)read_sz+/*padding*/(size_t)(read_sz%sizeof(MPI_Offset));

            nelems_to_read = (int)(read_sz/def_el_sz);
            if(var->dtype == MPI_FLOAT){
                double *ubuf;
                /*MPI_DOUBLE->MPI_FLOAT*/
                assert(ioreq->dtype == MPI_DOUBLE);
                ubuf = (double*)((unsigned char*)ioreq->user_buf+dbuf_offt);
                for(i = 0; i < nelems_to_read; i++){
                    ubuf[i] = (double)(((float*)cpbuf)[i]);
                    FARB_DBG(VERBOSE_ALL_LEVEL, "Recv: uval %.10f, bval %.10f",ubuf[i], ((float*)cpbuf)[i]);
                }
            } else if (var->dtype == MPI_DOUBLE){
                float *ubuf;
                /*MPI_FLOAT->MPI_DOUBLE*/
                assert(ioreq->dtype == MPI_FLOAT);
                ubuf = (float*)((unsigned char*)ioreq->user_buf+dbuf_offt);
                for(i = 0; i < nelems_to_read; i++){
                    ubuf[i] = (float)(((double*)cpbuf)[i]);
                    FARB_DBG(VERBOSE_ALL_LEVEL, "Recv: uval %.10f, bval %.10f",ubuf[i], ((double*)cpbuf)[i]);
                }
            } else {
                FARB_DBG(VERBOSE_ERROR_LEVEL, "This conversion type is not supported. Aborting.");
                assert(0);
            }
            farb_free(cpbuf, read_sz);

            ioreq->get_sz += (MPI_Offset)(read_sz/def_el_sz)*req_el_sz;
        }
        FARB_DBG(VERBOSE_ALL_LEVEL, "req %d, var %d, Got %d (expect %d)", ioreq->id, ioreq->var_id, (int)ioreq->get_sz, (int)ioreq->user_buf_sz);

        if(ioreq->get_sz == ioreq->user_buf_sz){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Complete req %d (left %d)", ioreq->id, fbuf->rreq_cnt-1);
            //delete this ioreq

            delete_ioreq(fbuf, &ioreq);


            if(fbuf->rreq_cnt == 0){
                FARB_DBG(VERBOSE_DBG_LEVEL, "Completed all rreqs. Notify master.");
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
    }
}

/*Send the file header and info about vars to the reader when the writer finishes the def mode*/
void send_file_info(file_buffer_t *fbuf, int reader_root)
{
    void *sbuf;
    MPI_Offset sbuf_sz, errno;

    FARB_DBG(VERBOSE_DBG_LEVEL, "Will send file info to reader %d", reader_root);
    pack_file_info(fbuf, &sbuf_sz, &sbuf);
    assert(sbuf_sz > 0);
    errno = MPI_Send(sbuf, (int)sbuf_sz, MPI_BYTE, reader_root, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm);
    CHECK_MPI(errno);
    gl_stats.nmsg_sent++;
    farb_free(sbuf, sbuf_sz);
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

    FARB_DBG(VERBOSE_DBG_LEVEL, "Mst %d will notify workgroup (msgtag %d) for %s", gl_my_rank, msgtag, fbuf->file_path);
    /*First, translate the ranks in the communicator which
    was used to open the file to global ranks*/
    MPI_Comm_group(gl_comps[gl_my_comp_id].comm, &glob_group);
    MPI_Comm_group(fbuf->comm, &file_group);

    nranks = fbuf->mst_info->my_workgroup_sz - 1;
    ranks = farb_malloc(nranks*sizeof(int));
    assert(ranks != NULL);
    glob_ranks = farb_malloc(nranks*sizeof(int));
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
    farb_free(ranks, nranks*sizeof(int));
    farb_free(glob_ranks, nranks*sizeof(int));
}

//TODO think about file overwriting
//file might not be closed
void progress_io_matching()
{
    MPI_Status status;
    int comp, flag, src, errno;
    file_buffer_t *fbuf;
    int bufsz;
    void *rbuf;
    char filename[MAX_FILE_NAME];
    int i, ncid;

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
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Recv root req from rank %d (comp %s)", src, gl_comps[comp].name);
                    rbuf = farb_malloc(MAX_FILE_NAME+2*sizeof(int));
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
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Don't know who is the root for %s now. Will queue the req.", fbuf->file_path);
                        file_info_req_q_t *req = farb_malloc(sizeof(file_info_req_q_t));
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
                        FARB_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");

                        memcpy(&fbuf->root_reader, (unsigned char*)rbuf+MAX_FILE_NAME, sizeof(int));
                        memcpy(&(fbuf->mst_info->nrranks_opened), (unsigned char*)rbuf+MAX_FILE_NAME+sizeof(int), sizeof(int));
                        send_file_info(fbuf, fbuf->root_reader);
                    } else {
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Forward the request to rank %d", fbuf->root_writer);
                        /*Forward this request to the root master*/
                        errno = MPI_Send(rbuf,(int)(MAX_FILE_NAME+2*sizeof(int)), MPI_BYTE, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[gl_my_comp_id].comm);
                        CHECK_MPI(errno);
                        gl_stats.nmatching_msg_sent++;
                        gl_stats.accum_msg_sz += (size_t)(MAX_FILE_NAME+2*sizeof(int));
                    }
                    farb_free(rbuf, MAX_FILE_NAME+2*sizeof(int));
                    break;
//                case FILE_READY_TAG:
//                    errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[comp].comm, &status);
//                    CHECK_MPI(errno);
//                    FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);
//
//                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//                    assert(fbuf != NULL);
//
//                    /*Receive the header*/
//                    MPI_Probe(src, HEADER_TAG, gl_comps[comp].comm, &status);
//                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
//                    fbuf->hdr_sz = (MPI_Offset)bufsz;
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Hdr size to receive %d", bufsz);
//                    fbuf->header = farb_malloc((size_t)bufsz);
//                    assert(fbuf->header != NULL);
//                    errno = MPI_Recv(fbuf->header, bufsz, MPI_BYTE, src, HEADER_TAG, gl_comps[comp].comm, &status);
//                    CHECK_MPI(errno);
//
//                    /*Receive info about vars*/
//                    MPI_Probe(src, VARS_TAG, gl_comps[comp].comm, &status);
//                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
//                    rbuf = farb_malloc(bufsz);
//                    assert(rbuf != NULL);
//                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, VARS_TAG, gl_comps[comp].comm, &status);
//                    CHECK_MPI(errno);
//                    unpack_vars(fbuf, bufsz, rbuf);
//                    farb_free(rbuf, bufsz);
//                    fbuf->is_ready = 1;
//                    break;
                case IO_REQS_TAG:
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Received reqs from %d, comp %d (my comp %d), bufsz %d", src, comp,gl_my_comp_id, (int)bufsz);
                    rbuf = farb_malloc(bufsz);
                    assert(rbuf != NULL);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_REQS_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    parse_ioreqs(rbuf, bufsz, src, gl_comps[comp].comm);
                    farb_free(rbuf, bufsz);
                    break;
//                case IO_REQ_TAG:
//                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Recved req src %d, comp %d (my comp %d), bufsz %d", status.MPI_SOURCE, comp,gl_my_comp_id, (int)bufsz);
//                    rbuf = farb_malloc(bufsz);
//                    assert(rbuf != NULL);
//                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_REQ_TAG, gl_comps[comp].comm, &status);
//                    CHECK_MPI(errno);
//                    parse_ioreq(rbuf, bufsz, src, gl_comps[comp].comm);
//                    farb_free(rbuf, bufsz);
//                    break;
                case IO_DATA_REQ_TAG:
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    rbuf = farb_malloc(bufsz);
                    assert(rbuf != NULL);

                    FARB_DBG(VERBOSE_DBG_LEVEL, "Recv data req from mst %d", src);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    send_data_wrt2rdr(rbuf, bufsz);
                    farb_free(rbuf, bufsz);
                    break;
                case IO_DATA_TAG:
                    //FARB_DBG(VERBOSE_DBG_LEVEL, "Recved data from %d", src);
                    MPI_Get_count(&status, MPI_BYTE, &bufsz);
                    rbuf = farb_malloc(bufsz);
                    assert(rbuf != NULL);
                    errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, IO_DATA_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);

                    recv_data_rdr(rbuf, bufsz);
                    farb_free(rbuf, bufsz);
                    break;
                case READ_DONE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, READ_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    assert(fbuf->root_writer == gl_my_rank);

                    if(comp == gl_my_comp_id){
                        fbuf->mst_info->nwranks_completed++;
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Recv read_done from writer %d (tot %d)", src, fbuf->mst_info->nwranks_completed);
                        /*shouldn't have matching going on for readers and writers at the same time*/
                        assert(fbuf->mst_info->nrranks_completed == 0);
                    } else {
                        assert(comp == fbuf->reader_id);
                        fbuf->mst_info->nrranks_completed++;
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Recv read_done from reader %d (tot %d)", src, fbuf->mst_info->nrranks_completed);
                        assert(fbuf->mst_info->nwranks_completed == 0);
                    }

                    if( ((fbuf->mst_info->nrranks_opened > 0) && (fbuf->mst_info->nrranks_completed == fbuf->mst_info->nrranks_opened)) ||
                        ((fbuf->mst_info->nwranks_opened > 0) && (fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened)) ){

                            FARB_DBG(VERBOSE_DBG_LEVEL, "Notify masters that matching for %s is finished", fbuf->file_path);
                            for(i = 0; i < fbuf->mst_info->nmasters; i++){
                                    if(fbuf->mst_info->masters[i] == gl_my_rank)
                                        continue;
                                    FARB_DBG(VERBOSE_DBG_LEVEL, "Notify mst %d", fbuf->mst_info->masters[i]);
                                    errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], MATCH_DONE_TAG, gl_comps[gl_my_comp_id].comm);
                                    CHECK_MPI(errno);
                                    gl_stats.nmsg_sent++;
                            }

                            if(fbuf->mst_info->is_master_flag){
                                FARB_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                                /*Tell other writer ranks that they can complete matching*/
                                notify_workgroup(fbuf, MATCH_DONE_TAG);
                                FARB_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                                fbuf->done_matching_flag = 1;
                                if(fbuf->rdr_closed_flag){
                                    FARB_DBG(VERBOSE_DBG_LEVEL, "Cleaning up all reqs and dbs");
                                    assert(fbuf->rreq_cnt == 0);
                                    /*Complete my own write requests*/
                                    delete_ioreqs(fbuf);

                                    if(fbuf->mst_info->is_master_flag){
                                        /*Check that I don't have any read reqs incompleted*/
                                        assert(fbuf->mst_info->iodb->nritems == 0);
                                        /*Clean my iodb*/
                                        clean_iodb(fbuf->mst_info->iodb);
                                    }
                                }
                            }
                    }

//
//                    if(fbuf->mst_info->is_master_flag){ //I am master
//
//                        if(comp == gl_my_comp_id){
//                            if(src == fbuf->root_writer){
//                                FARB_DBG(VERBOSE_DBG_LEVEL, "Recevied DONE from mst 0");
//                                //set value to notify all writers
//                                fbuf->mst_info->nwranks_completed = fbuf->mst_info->nwranks_opened;
//                            } else {
//                                assert(gl_my_rank == fbuf->root_writer);
//                                fbuf->mst_info->nwranks_completed++;
//                                FARB_DBG(VERBOSE_DBG_LEVEL, "Recevied DONE from writer %d (tot done %u)", status.MPI_SOURCE, fbuf->mst_info->nwranks_completed);
//
//                                if(fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened){
//                                /*Notify other masters that all writers finished reading the file*/
//                                  FARB_DBG(VERBOSE_DBG_LEVEL, "Notify masters that writers finished reading");
//                                  for(i = 0; i < fbuf->mst_info->nmasters; i++){
//                                    if(fbuf->mst_info->masters[i] == gl_my_rank)
//                                        continue;
//                                    errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], MATCH_DONE_TAG, gl_comps[gl_my_comp_id].comm);
//                                    CHECK_MPI(errno);
//                                  }
//                                }
//                            }
//                         }else{
//                            //Only master 0 is supposed to receive this message from a reader rank
//                            assert(gl_my_rank == fbuf->root_writer);
//                            fbuf->mst_info->nrranks_completed++; //this message was from a reader
//                            FARB_DBG(VERBOSE_DBG_LEVEL, "Recevied DONE from reader %d (tot done %u)", status.MPI_SOURCE, fbuf->mst_info->nrranks_completed);
//
//                            if(fbuf->mst_info->nrranks_completed == fbuf->mst_info->nrranks_opened){
//                                /*Notify other masters that readers have finished */
//                                  FARB_DBG(VERBOSE_DBG_LEVEL, "Notify masters readers completed");
//                                  for(i = 0; i < fbuf->mst_info->nmasters; i++){
//                                    if(fbuf->mst_info->masters[i] == gl_my_rank)
//                                        continue;
//                                    errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], MATCH_DONE_TAG, gl_comps[gl_my_comp_id].comm);
//                                    CHECK_MPI(errno);
//                                  }
//                            }
//                         }


//                        if( ((fbuf->mst_info->nrranks_opened > 0) && (fbuf->mst_info->nrranks_completed == fbuf->mst_info->nrranks_opened)) ||
//                            ((fbuf->mst_info->nwranks_opened > 0) && (fbuf->mst_info->nwranks_completed == fbuf->mst_info->nwranks_opened)) ){
//
//                            FARB_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
//                            /*Tell other writer ranks that they can complete matching*/
//                            notify_workgroup(fbuf, MATCH_DONE_TAG);
//                            FARB_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
//                            fbuf->done_matching_flag = 1;
//                            if(fbuf->fclosed_flag){
//                                /*Check that I don't have any read reqs incompleted*/
//                                assert(fbuf->mst_info->iodb->nritems == 0);
//                                assert(fbuf->rreq_cnt == 0);
//                                /*Complete my own write requests*/
//                                delete_ioreqs(fbuf);
//                                /*Clean my iodb*/
//                                clean_iodb(fbuf->mst_info->iodb);
//                            }
//                        }
//
//                    } else { /*I am writer*/
//                        FARB_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
//                        fbuf->done_matching_flag = 1;
//                        if(fbuf->fclosed_flag){
//                            assert(fbuf->rreq_cnt == 0);
//                            /*It's safe to delete all wreqs now.*/
//                            delete_ioreqs(fbuf);
//                        }
//                    }

                    break;
                case MATCH_DONE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, MATCH_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    if(fbuf->mst_info->is_master_flag){
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                        /*Tell other writer ranks that they can complete matching*/
                        notify_workgroup(fbuf, MATCH_DONE_TAG);
                    }

                    FARB_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                    fbuf->done_matching_flag = 1;
                    if(fbuf->rdr_closed_flag){
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Cleaning up all reqs and dbs");
                        assert(fbuf->rreq_cnt == 0);
                        /*Complete my own write requests*/
                        delete_ioreqs(fbuf);

                        if(fbuf->mst_info->is_master_flag){
                            /*Check that I don't have any read reqs incompleted*/
                            assert(fbuf->mst_info->iodb->nritems == 0);
                            /*Clean my iodb*/
                            clean_iodb(fbuf->mst_info->iodb);
                        }
                    }
                    break;
                case IO_CLOSE_FILE_TAG:
                    errno = MPI_Recv(&ncid, 1, MPI_INT, src, IO_CLOSE_FILE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
                    assert(fbuf != NULL);
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Recv close file tag for %s from %d", fbuf->file_path, src);
                    assert(!fbuf->rdr_closed_flag); //check that haven't received that before

                    if(gl_my_rank == fbuf->root_writer){
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Notify other masters that readers are closing the file");
                        for(i = 0; i < fbuf->mst_info->nmasters; i++){
                            if(fbuf->mst_info->masters[i] == gl_my_rank)
                                continue;
                            errno = MPI_Send(&ncid, 1, MPI_INT, fbuf->mst_info->masters[i], IO_CLOSE_FILE_TAG, gl_comps[gl_my_comp_id].comm);
                            CHECK_MPI(errno);
                            gl_stats.nmsg_sent++;
                        }
                    }
                    if(fbuf->mst_info->is_master_flag){
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Notify writers that they can close the file %s", fbuf->file_path);
                        notify_workgroup(fbuf, IO_CLOSE_FILE_TAG);
                    }

                    fbuf->rdr_closed_flag = 1;
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Close flag set for file %s", fbuf->file_path);
                    if(fbuf->done_matching_flag){
                        if(fbuf->mst_info->is_master_flag){
                            /*Check that I don't have any read reqs incompleted*/
                            assert(fbuf->mst_info->iodb->nritems == 0);
                            /*Clean my iodb*/
                            clean_iodb(fbuf->mst_info->iodb);
                        }
                         /*Can delete write requests only if all ps-s have finished
                         matching.*/
                        delete_ioreqs(fbuf);
                    }
                    break;
                case ROOT_MST_TAG:
                    src = status.MPI_SOURCE;
                    errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, ROOT_MST_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(errno);
                    FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive ROOT_MST_TAG notif for %s", filename);

                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    assert(fbuf->root_writer == -1);
                    fbuf->root_writer = src;
                    break;
                default:
                    FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: unknown tag %d", status.MPI_TAG);
                    assert(0);
            }
        }
    }
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
//            if(fbuf->mst_info->nmasters == 0)
//                fbuf->mst_info->nmasters++;
//            else if(nranks % wg > (int)(wg/2) ){
//                if(myrank >= fbuf->mst_info->nmasters * wg){
//                    fbuf->mst_info->my_master = fbuf->mst_info->nmasters * wg;
//                    gl_conf.my_workgroup_sz = nranks % wg;
//                }
//                fbuf->mst_info->nmasters++;
//            } else if ( (nranks % wg > 0) && (myrank >= (fbuf->mst_info->nmasters-1)*wg)){
//                /*Merge last smaller group with the previous group*/
//                fbuf->mst_info->my_master = (fbuf->mst_info->nmasters-1) * wg;
//                gl_conf.my_workgroup_sz = wg + nranks % wg;
//            }
    }

    if(myrank == 0)
        FARB_DBG(VERBOSE_ALL_LEVEL, "Nmasters %d", mst_info->nmasters);

    mst_info->masters = (int*)farb_malloc(mst_info->nmasters * sizeof(int));
    assert(mst_info->masters != NULL);
    masters = (int*)farb_malloc(mst_info->nmasters * sizeof(int));
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

//    if(gl_my_rank == my_master_glob)
//        mst_info->is_master_flag = 1;
//    else
//        mst_info->is_master_flag = 0;
    mst_info->is_master_flag = (gl_my_rank == my_master_glob) ? 1 : 0;

    if(myrank == 0){
        for(i = 0; i < mst_info->nmasters; i++)
            FARB_DBG(VERBOSE_DBG_LEVEL, "Rank %d is a master", mst_info->masters[i]);
    }
    if(mst_info->is_master_flag)
        FARB_DBG(VERBOSE_DBG_LEVEL, "My wg size %d", mst_info->my_workgroup_sz);
    farb_free(masters, mst_info->nmasters * sizeof(int));
    MPI_Group_free(&glob_group);
    MPI_Group_free(&file_group);


//    //For each component, find out to which master I should send read requests
//    for(i=0; i < gl_ncomp; i++){
//        if(i == gl_my_comp_id)
//            continue;
//        MPI_Comm_remote_size(gl_comps[i].comm, &nrranks);
//
//        if(myrank < nrranks)
//            gl_comps[i].master = fbuf->mst_info->my_master; //use same rank
//        else {
//            int nmasters = (int)(nrranks/wg);
//            if( nrranks % wg > (int)(wg/2))
//               nmasters++;
//            gl_comps[i].master = (nmasters-1)*wg;
//        }
//    }

    return 0;
}

