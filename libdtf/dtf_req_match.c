#include "dtf_req_match.h"
#include "dtf_util.h"
#include "dtf.h"
#include <unistd.h>


static void send_data_wrt2rdr(void* buf, int bufsz);

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


void delete_ioreq(file_buffer_t *fbuf, io_req_t **ioreq)
{
	DTF_DBG(VERBOSE_ALL_LEVEL, "Delete req %d, cur wreqs %d, rreqs %d", (*ioreq)->id,
			fbuf->wreq_cnt, fbuf->rreq_cnt);

    dtf_var_t *var = fbuf->vars[(*ioreq)->var_id];
    
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
    DTF_DBG(VERBOSE_ALL_LEVEL, "Deleted, cur wreqs %d, rreqs %d",
			fbuf->wreq_cnt, fbuf->rreq_cnt);
}

void delete_ioreqs(file_buffer_t *fbuf, int finalize)
{
    io_req_t *ioreq, *tmp;
    DTF_DBG(VERBOSE_ALL_LEVEL, "Delete io requests for file %s", fbuf->file_path);

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
		if(ioreq->is_permanent && !finalize){
			ioreq = ioreq->next;
			continue;
		}
		tmp = ioreq->next;
        delete_ioreq(fbuf, &ioreq);
        ioreq = tmp; //fbuf->ioreqs;
    }

    if(finalize){
		assert(fbuf->rreq_cnt == 0);
		assert(fbuf->wreq_cnt == 0);
		assert(fbuf->ioreqs == NULL);
	}
}

static void pack_file_info(file_buffer_t *fbuf, MPI_Offset *bufsz, void **buf)
{
    dtf_var_t *var;
    MPI_Offset sz = 0, offt = 0;
    unsigned char *chbuf;
    int i;
 //   rb_red_blk_node *var_node;
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

    for(i = 0; i < fbuf->nvars; i++){
//        tmp_var.id = i;
//        var_node = RBExactQuery(fbuf->vars, &tmp_var);
//        assert(var_node != NULL);
//        var = (dtf_var_t*)(var_node->key);
        /*sz += sizeof(var->id) + sizeof(var->el_sz) +
        sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims;*/
        sz += sizeof(MPI_Offset)*3+ sizeof(MPI_Offset)*fbuf->vars[i]->ndims;
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
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->nvars;
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d vars", fbuf->nvars);
    /*vars*/

    for(i = 0; i < fbuf->nvars; i++){
        var = fbuf->vars[i];
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
    int i, varid, nvars;
    file_buffer_t *fbuf;
    int ndims;
    dtf_var_t *var;
    MPI_Offset offt = 0;
    int type;
    MPI_Datatype dtype;
    MPI_Offset *shape;
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
    nvars = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", nvars);
    /*vars*/
    for(i = 0; i < nvars; i++){
        varid = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        assert(varid >= fbuf->nvars);

        type = (int)(*((MPI_Offset*)(chbuf+offt)));
        dtype = int2mpitype(type);
        offt+=sizeof(MPI_Offset);
        ndims = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        shape = (MPI_Offset*)(chbuf+offt);
        offt += sizeof(MPI_Offset)*ndims;

        var = new_var(varid, ndims, dtype, shape);
        add_var(fbuf, var);
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
    block_t *wblock;
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

    double t_st, t_idle = 0, t_start;

    int n_matched_blocks = 0;
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
    nwriters = 0;
    ritem = fbuf->mst_info->iodb->ritems;
    while(ritem != NULL){
		
        t_st = MPI_Wtime();

        /*If we are supposed to match intra component io requests
        then skip all other requests*/
        if(intracomp_io_flag && (ritem->comm != gl_comps[gl_my_comp_id].comm)){
            gl_stats.idle_time += MPI_Wtime() - t_st;
            t_idle += MPI_Wtime() - t_st;
            gl_stats.idle_do_match_time += MPI_Wtime() - t_st;
            ritem = ritem->next;
            continue;
        }

        n_matched_blocks = 0;

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
            t_st = MPI_Wtime();
            
            witem = fbuf->mst_info->iodb->witems[var_id];

			if(witem == NULL){
				rblock = rblock->next;
				/*No match right now*/
				gl_stats.idle_time += MPI_Wtime() - t_st;
				t_idle += MPI_Wtime() - t_st;
				gl_stats.idle_do_match_time += MPI_Wtime() - t_st;
				continue;
			}
                       
            var = fbuf->vars[var_id];          
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
				if(var->ndims > 0)
					wblock = rb_find_block(witem->dblocks, rblock->start, var->ndims);
				else
					wblock = (block_t*)witem->dblocks;

                if(wblock == NULL){

                    //didn't find
                    gl_stats.idle_time += MPI_Wtime() - t_st;
                    t_idle += MPI_Wtime() - t_st;
                    gl_stats.idle_do_match_time += MPI_Wtime() - t_st;
                    break;
                }

				DTF_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for var %d agains block:", var->id);
				for(i = 0; i < var->ndims; i++)
					DTF_DBG(VERBOSE_DBG_LEVEL, "%lld -> %lld || %lld ->%lld",
							rblock->start[i], rblock->count[i], wblock->start[i], wblock->count[i]);
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
                t_st = MPI_Wtime(); //reset
            }

            dtf_free(matched_count, ndims*sizeof(MPI_Offset));

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
    if(nwriters == 0)
        t_st = MPI_Wtime(); //reset
    /*dealloc stuff*/
    dtf_free(writers, allocd_nwriters*sizeof(int));
    dtf_free(offt, allocd_nwriters*sizeof(size_t));
    dtf_free(bufsz, allocd_nwriters*sizeof(int));
    dtf_free(sbuf, allocd_nwriters*sizeof(unsigned char*));
    gl_stats.ndb_match++;

    if(nwriters == 0)
        gl_stats.idle_time += MPI_Wtime() - t_st;

    assert( MPI_Wtime() - t_start >= t_idle);
    gl_stats.master_time = MPI_Wtime() - t_start - t_idle;  //useful work


    DTF_DBG(VERBOSE_DBG_LEVEL, "after matching: %d ritems", (int)fbuf->mst_info->iodb->nritems);
}

void init_iodb(file_buffer_t *fbuf)
{
    fbuf->mst_info->iodb = dtf_malloc(sizeof(struct ioreq_db));
    assert(fbuf->mst_info->iodb != NULL);
    fbuf->mst_info->iodb->ritems = NULL;
    fbuf->mst_info->iodb->witems = NULL;   
    fbuf->mst_info->iodb->nritems = 0;
    fbuf->mst_info->iodb->updated_flag = 0;
}


void clean_iodb(ioreq_db_t *iodb, int nvars)
{
	int i;
	write_db_item_t *witem;
	read_db_item_t *ritem;
	
    if(iodb == NULL){
		return;
	}
//    DTF_DBG(VERBOSE_DBG_LEVEL, "Clean iodb: memuse %lu, peak %lu", gl_stats.iodb_cur_memuse, gl_stats.iodb_peak_memuse);

	if(iodb->witems != NULL){
		for(i = 0; i < nvars; i++){
		
			if(iodb->witems[i] == NULL)
				continue;
			
			witem = iodb->witems[i];
			
			if(witem->ndims > 0){
				RBTreeDestroy(witem->dblocks);
				gl_stats.malloc_size -= witem->nblocks*(sizeof(block_t)+sizeof(MPI_Offset)*2*witem->ndims);
			} else
				dtf_free(witem->dblocks, sizeof(block_t));
				
			dtf_free(witem, sizeof(write_db_item_t));
			iodb->witems[i] = NULL;
		}
	}

	ritem = iodb->ritems;
    while(ritem != NULL){
        if(ritem->dblocks != NULL){
            read_dblock_t *block = ritem->dblocks;
            while(block != NULL){

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
  
    iodb->updated_flag = 0;
   
	DTF_DBG(VERBOSE_DBG_LEVEL, "iodb clean");
}


void finalize_iodb(file_buffer_t *fbuf)
{
	if(fbuf->mst_info->iodb == NULL)
		return;
	
	clean_iodb(fbuf->mst_info->iodb, fbuf->nvars);
	dtf_free(fbuf->mst_info->iodb->witems, fbuf->nvars*sizeof(write_db_item_t*));
	dtf_free(fbuf->mst_info->iodb, sizeof(ioreq_db_t));
	fbuf->mst_info->iodb = NULL;
}


static void parse_ioreqs(void *buf, int bufsz, int rank, MPI_Comm comm)
{
    int var_id, rw_flag;
    dtf_var_t *var = NULL;
    file_buffer_t *fbuf;
    size_t offt = 0;
    char filename[MAX_FILE_NAME];
    unsigned char *chbuf = (unsigned char*)buf;
	read_db_item_t *ritem = NULL;
    memcpy(filename, chbuf, MAX_FILE_NAME);
    offt += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    assert(fbuf != NULL);
    
    double t_st = MPI_Wtime();
    
    if( (gl_conf.iodb_build_mode == IODB_BUILD_RANK) && fbuf->done_matching_flag){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: discard new ioreqs from %d because already finished matching", rank);
		return;
	}
    
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

        gl_stats.iodb_nioreqs++;
		var = fbuf->vars[var_id];
		
        if(rw_flag == DTF_READ){
			 
			 if(ritem == NULL){
			
				ritem = fbuf->mst_info->iodb->ritems;
				while(ritem != NULL){
					if( (ritem->rank == rank) && (ritem->comm == comm) )
						break;
					ritem = ritem->next;
				}

				if(ritem == NULL){
					ritem = (read_db_item_t*)dtf_malloc(sizeof(read_db_item_t));
					assert(ritem != NULL);
					ritem->rank = rank;
					ritem->comm = comm;
					ritem->dblocks = NULL;
					ritem->nblocks = 0;
					ritem->last_block = NULL;
					ritem->next = NULL;
					ritem->prev = NULL;
					//enqueue
					if(fbuf->mst_info->iodb->ritems == NULL)
						fbuf->mst_info->iodb->ritems = ritem;
					else{
						ritem->next = fbuf->mst_info->iodb->ritems;
						fbuf->mst_info->iodb->ritems->prev = ritem;
						fbuf->mst_info->iodb->ritems = ritem;
					}
					fbuf->mst_info->iodb->nritems++;
				} 
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
            if(ritem->dblocks == NULL){
                ritem->dblocks = dblock;
                ritem->last_block = dblock;
            } else {
                ritem->last_block->next = dblock;
                dblock->prev = ritem->last_block;
                ritem->last_block = dblock;
            }
            ritem->nblocks++;
            DTF_DBG(VERBOSE_DBG_LEVEL, "ritems %d, cur item (r %d) %lld blocks",(int)fbuf->mst_info->iodb->nritems, ritem->rank, ritem->nblocks);
        } else { /*DTF_WRITE*/
            write_db_item_t *witem;
            /*Allow write requests only from the writer*/
            assert(comm == gl_comps[fbuf->writer_id].comm);

            if(fbuf->mst_info->iodb->witems[var_id] == NULL){
                fbuf->mst_info->iodb->witems[var_id] = (write_db_item_t*)dtf_malloc(sizeof(write_db_item_t));
                assert(fbuf->mst_info->iodb->witems[var_id] != NULL);
                witem = fbuf->mst_info->iodb->witems[var_id];
               
                if(var->ndims > 0)
                    witem->dblocks = (void*)RBTreeCreateBlocks(rb_key_cmp, NullFunction, rb_destroy_node_info, rb_print_key, rb_print_info, 0);
                else{
                    witem->dblocks = dtf_malloc(sizeof(block_t));
                    assert(witem->dblocks != NULL);
				}
                witem->nblocks = 0;
                witem->ndims = var->ndims;
            } else 
				witem = fbuf->mst_info->iodb->witems[var_id];
            
            if(var->ndims > 0){
				insert_info *info = dtf_malloc(sizeof(insert_info));
				assert(info != NULL);
				info->ndims = var->ndims;

				info->blck = dtf_malloc(sizeof(block_t));
				assert(info->blck != NULL);
				info->cur_dim = 0;
				info->blck->rank = rank;
				info->blck->start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                assert(info->blck->start != NULL);
                info->blck->count = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                assert(info->blck->count != NULL);
                memcpy(info->blck->start, chbuf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);
                memcpy(info->blck->count, chbuf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);
                /*add block to the database*/
				rb_red_blk_node *bl_node = RBTreeInsertBlock(witem->dblocks, info);
				assert(bl_node != NULL);
				dtf_free(info, sizeof(insert_info));

			} else{
                ((block_t*)witem->dblocks)->start = NULL;
                ((block_t*)witem->dblocks)->count = NULL;
                ((block_t*)witem->dblocks)->rank = rank;
            }
            witem->nblocks++;
        }
    }
    assert(offt == (size_t)bufsz);
	fbuf->mst_info->iodb->updated_flag = 1;
    gl_stats.parse_ioreq_time += MPI_Wtime() - t_st;
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
        char *s = getenv("BT_BENCH");
        int bt_bench = 0;
        if(s!=NULL)
            bt_bench = 1;
        if(bt_bench){
            int def_el_sz;
            MPI_Type_size(MPI_DOUBLE, &def_el_sz);
            ioreq->user_buf_sz *= def_el_sz;
        } else
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
    ioreq->is_permanent = 0;

    if( (rw_flag == DTF_WRITE) && gl_conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
        ioreq->checksum = compute_checksum(buf, ndims, count, dtype);
        DTF_DBG(VERBOSE_DBG_LEVEL, "checksum %.4f", ioreq->checksum);
    } else
        ioreq->checksum = 0;
    gl_stats.nioreqs++;
    return ioreq;
}

//TODO for writing scalar var only one ps will save an ioreq
/*The metadata is distributed among masters based on the var id.*/
void send_ioreqs_by_var(file_buffer_t *fbuf, int intracomp_match)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int mst = 0;
    int nrreqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    int data_to_send = 0;
    int nmasters = fbuf->mst_info->nmasters;
    int idx, err, flag;
    MPI_Request *reqs;
    
    if(fbuf->ioreqs == NULL)
        return;
    reqs = dtf_malloc(nmasters * sizeof(MPI_Request));
    assert(reqs != NULL);
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
            //All the following reqs should have been sent already
            break;
        }
        data_to_send = 1;

        if( (fbuf->writer_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ) && intracomp_match)
            nrreqs++;
        else if((fbuf->reader_id == gl_my_comp_id) && (ioreq->rw_flag == DTF_READ))
            nrreqs++;

        var = fbuf->vars[ioreq->var_id];
        mst = ioreq->var_id % fbuf->mst_info->nmasters;
        bufsz[mst] += sizeof(MPI_Offset)*2 + var->ndims*2*sizeof(MPI_Offset);
        ioreq = ioreq->next;
    }
    if(nrreqs == 0){
        if((gl_my_comp_id == fbuf->writer_id) && intracomp_match){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
            dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
            int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer,
                      READ_DONE_TAG, gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_msg_q);
        } else if(gl_my_comp_id == fbuf->reader_id){
            if(gl_my_rank == fbuf->root_reader){
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
				err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_msg_q);
			} else {
				err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm);
				CHECK_MPI(err);
				fbuf->done_matching_flag = 1;
			}
		}
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
            //All the following reqs should have been sent already
            break;
        }
        var = fbuf->vars[ioreq->var_id];
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
    flag = 0;
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


/*This version divides the variable data among masters along the var's biggest dimension.
  If there is an unlimited dimension, the division happens along it*/
void send_ioreqs_by_block(file_buffer_t *fbuf, int intracomp_match)
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
        } else if(gl_my_comp_id == fbuf->reader_id){
            if(gl_my_rank == fbuf->root_reader){
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
				err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_msg_q);
			} else {
				err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm);
				CHECK_MPI(err);
				fbuf->done_matching_flag = 1;
			}
		}
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            break;
        }
        var = fbuf->vars[ioreq->var_id];
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
//            if(gl_conf.block_sz_range == AUTO_BLOCK_SZ_RANGE){ //divide equally among masters by the 1st dimension
//                //If the first dimension is unlimited then will divide by default blocks
//                if(var->shape[0] == DTF_UNLIMITED)
//                    block_range = DEFAULT_BLOCK_SZ_RANGE;
//                else{
//                    block_range = (MPI_Offset)(var->shape[0]/fbuf->mst_info->nmasters);
//                    if(block_range == 0)
//                        block_range = var->shape[0];
//                }
//            } else
//                block_range = gl_conf.block_sz_range;
            if(var->shape[var->max_dim] == DTF_UNLIMITED)
                block_range = DEFAULT_BLOCK_SZ_RANGE;
            else{
                block_range = (MPI_Offset)(var->shape[var->max_dim]/fbuf->mst_info->nmasters);
                if(var->shape[var->max_dim]%fbuf->mst_info->nmasters>0)
					block_range++;
			}
            if(block_range == 0)
                block_range = DEFAULT_BLOCK_SZ_RANGE;

            DTF_DBG(VERBOSE_DBG_LEVEL, "Var %d, max dim %d, block_range %lld", var->id, var->max_dim, block_range);

            shift = 0;
            while(ioreq->start[var->max_dim] + shift < ioreq->start[var->max_dim] + ioreq->count[var->max_dim]){

                if(block_range == DEFAULT_BLOCK_SZ_RANGE)
                    mst = (int)((ioreq->start[var->max_dim] + shift) % fbuf->mst_info->nmasters);
                else
                    mst = (int)((ioreq->start[var->max_dim] + shift)/block_range);
                DTF_DBG(VERBOSE_DBG_LEVEL, "mst %d", mst);
                assert(mst < fbuf->mst_info->nmasters);
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

                /*Adjust corresponding start coordinate*/
                *(MPI_Offset*)(sbuf[mst]+offt[mst] + var->max_dim*sizeof(MPI_Offset)) = ioreq->start[var->max_dim] + shift;
                strt = (MPI_Offset*)(sbuf[mst]+offt[mst]);
                offt[mst] += var->ndims*sizeof(MPI_Offset);
                memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));

                /*Adjust corresponding count*/
                if(ioreq->count[var->max_dim] - shift > block_range )
                    *(MPI_Offset*)(sbuf[mst]+offt[mst] + var->max_dim*sizeof(MPI_Offset)) = block_range;
                else
                    *(MPI_Offset*)(sbuf[mst]+offt[mst]+var->max_dim*sizeof(MPI_Offset)) = ioreq->count[var->max_dim] - shift;
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
            int flag;
            MPI_Status status;
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], IO_REQS_TAG);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu)", fbuf->mst_info->masters[mst], offt[mst], bufsz[mst]);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, fbuf->mst_info->masters[mst], IO_REQS_TAG,
                            gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            err = MPI_Test(&(msg->req), &flag, &status);
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

/*Writers send their ioreqs to every master. Reader sends all its ioreqs to 
 * its designated master*/
void send_ioreqs_by_mst(file_buffer_t *fbuf, int intracomp_match)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int err;
    int nrreqs = 0;
    unsigned char *sbuf = NULL;
    size_t bufsz = 0, offt = 0;
    int data_to_send = 0;

    if(fbuf->ioreqs == NULL)
        return;   

	//count bufsz to allocate
    bufsz += MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
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

        var = fbuf->vars[ioreq->var_id];
        
        bufsz += sizeof(MPI_Offset)*2 + var->ndims*2*sizeof(MPI_Offset);
        ioreq = ioreq->next;
    }
    
    if(nrreqs == 0){
        if((gl_my_comp_id == fbuf->writer_id) && intracomp_match){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
            dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
            int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer,
                      READ_DONE_TAG, gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_msg_q);
        } else if(gl_my_comp_id == fbuf->reader_id){
            if(gl_my_rank == fbuf->root_reader){
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
				err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_msg_q);
			} else {
				err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm);
				CHECK_MPI(err);
				fbuf->done_matching_flag = 1;
			}
		}
    }
    if(!data_to_send)
        return;   //nothing to send

	DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu", bufsz);
	sbuf = dtf_malloc(bufsz);
	assert(sbuf != NULL);

	memcpy(sbuf, fbuf->file_path, MAX_FILE_NAME);
	offt = MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
    
    //pack ioreqs
    ioreq = fbuf->ioreqs;
    while(ioreq != NULL){
        if(ioreq->sent_flag){
            //All the following reqs should have been sent already
            break;
        }
        var = fbuf->vars[ioreq->var_id];
       
        /*Store var_id, rw_flag, start[] and count[]*/
        *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->rw_flag;
        offt += sizeof(MPI_Offset);
        *(MPI_Offset*)(sbuf + offt) = (MPI_Offset)ioreq->var_id;
        offt += sizeof(MPI_Offset);
        memcpy(sbuf+offt, ioreq->start, var->ndims*sizeof(MPI_Offset));
        offt += var->ndims*sizeof(MPI_Offset);
        memcpy(sbuf+offt, ioreq->count, var->ndims*sizeof(MPI_Offset));
        offt += var->ndims*sizeof(MPI_Offset);
        ioreq->sent_flag = 1;
        ioreq = ioreq->next;
    }
    
	assert(offt == bufsz);

    //send ioreqs
	if(gl_my_comp_id == fbuf->writer_id){
		int mine = 0, mst;
				
		for(mst = 0; mst < fbuf->mst_info->nmasters; mst++){
			
			if( fbuf->mst_info->masters[mst] == gl_my_rank){
				mine = 1;
			} else {
				int flag;
				MPI_Status status;
				dtf_msg_t *msg;
				void *buf;
				
				buf = dtf_malloc(bufsz);
				assert(buf != NULL);
			    memcpy(buf, sbuf, bufsz);
			    
				msg = new_dtf_msg(buf, bufsz, IO_REQS_TAG);
				DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu)", fbuf->mst_info->masters[mst], offt, bufsz);
				err = MPI_Isend(buf, (int)offt, MPI_BYTE, fbuf->mst_info->masters[mst], IO_REQS_TAG,
								gl_comps[fbuf->writer_id].comm, &(msg->req));
				CHECK_MPI(err);
				err = MPI_Test(&(msg->req), &flag, &status);
				CHECK_MPI(err);
				ENQUEUE_ITEM(msg, gl_msg_q);
			}
		}
		if(mine)
			parse_ioreqs(sbuf, (int)offt, gl_my_rank, gl_comps[gl_my_comp_id].comm);
		
	} else {
		int my_mst = gl_my_rank % fbuf->mst_info->nmasters;
        assert(gl_my_comp_id == fbuf->reader_id);
        DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu)", fbuf->mst_info->masters[my_mst], offt);
		err = MPI_Send((void*)sbuf, (int)offt, MPI_BYTE, fbuf->mst_info->masters[my_mst], IO_REQS_TAG,
							gl_comps[fbuf->writer_id].comm);
		CHECK_MPI(err);
	
	}
	dtf_free(sbuf, bufsz);
}



static void reset_fbuf_match(file_buffer_t *fbuf)
{
    if ((fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag)
        /*Reset flag for future matchings*/
        fbuf->mst_info->iodb->updated_flag = 0;

    fbuf->is_matching_flag = 0;
    fbuf->mst_info->nrranks_completed = 0;
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


void skip_match(file_buffer_t *fbuf, const char *filename, MPI_Comm comm, int writer_id)
{
	int err;
	
	if( fbuf == NULL){
		int rank;
		
		MPI_Comm_rank(comm, &rank);
        if(rank == 0){
			char fpath[MAX_FILE_NAME];
			int root_writer = inquire_root(filename);
			
			DTF_DBG(VERBOSE_DBG_LEVEL, "Send skip notif for file %s to writer %d", filename, root_writer);
			err = MPI_Send((void*)filename, MAX_FILE_NAME, MPI_CHAR, root_writer, SKIP_MATCH_TAG, gl_comps[writer_id].comm);
			CHECK_MPI(err);
			DTF_DBG(VERBOSE_DBG_LEVEL, "Wait for confirmation");
			err = MPI_Recv(fpath, MAX_FILE_NAME, MPI_CHAR, root_writer, READ_DONE_CONFIRM_TAG, gl_comps[writer_id].comm, MPI_STATUS_IGNORE);
			CHECK_MPI(err);
			assert(strcmp(fpath, filename)==0);
			
		}
		//Call barrier to make sure that we do not mix up dtf_skip_match with 
		//any further dtf_match_io calls
		MPI_Barrier(comm);
       
	} else {
		//Treat this as a normal match I/O where reader doesn't have any read ioreqs
		
		//First check that received confirmation from the writer
		// from any previous matches
		if(fbuf->done_match_confirm_flag == DTF_UNDEFINED)
			fbuf->done_match_confirm_flag = 0;
		else{
			if(gl_my_rank == fbuf->root_reader)
				while(!fbuf->done_match_confirm_flag){
					progress_io_matching();
				}

			int err = MPI_Bcast(&(fbuf->done_match_confirm_flag), 1, MPI_INT, 0, fbuf->comm);
			CHECK_MPI(err);
			assert(fbuf->done_match_confirm_flag);
			
			fbuf->done_match_confirm_flag = 0; //reset
		}
		
		if(gl_my_rank == fbuf->root_reader){
			int i;
			
			DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer masters readers completed (skip) matching");
			for(i = 0; i < fbuf->mst_info->nmasters; i++){
				int err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->masters[i], READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
				CHECK_MPI(err);
			}
		}
	}
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Skip match done");
}

int match_ioreqs(file_buffer_t *fbuf, int intracomp_io_flag)
{

//    double t_all = MPI_Wtime();
//    double t_part1 = MPI_Wtime();
//    double t_part2=MPI_Wtime();
//    double t_part3=MPI_Wtime();
    double t_start;
    double t_st;
    int finalize = 0;
	t_start = MPI_Wtime();
    DTF_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs for file %s (ncid %d), intracomp %d", fbuf->file_path, fbuf->ncid, intracomp_io_flag);

    if(fbuf->mst_info->is_master_flag){
        //there should be no unfinished matches
        assert(fbuf->mst_info->nrranks_completed == 0);
    }
    assert(!fbuf->is_matching_flag);

    if(gl_my_comp_id == fbuf->reader_id){
		/*Reader cannot proceed to a new match for this file until
		 * writer confirmed that it has finished with the previous match.
		 * Crucial for multi-iterative programs where there is I/O in
		 * every iteration*/
		if(fbuf->done_match_confirm_flag == DTF_UNDEFINED)
			fbuf->done_match_confirm_flag = 0;
		else{
			if(gl_my_rank == fbuf->root_reader)
				while(!fbuf->done_match_confirm_flag){
					progress_io_matching();
				}

			int err = MPI_Bcast(&(fbuf->done_match_confirm_flag), 1, MPI_INT, 0, fbuf->comm);
			CHECK_MPI(err);
			assert(fbuf->done_match_confirm_flag);
			DTF_DBG(VERBOSE_DBG_LEVEL, "Confirmation received. Proceed to match");
			fbuf->done_match_confirm_flag = 0; //reset
		}
	} else 
		progress_io_matching();

    fbuf->done_matching_flag = 0;
    
    /*If a writer process doesn't have any io requests, it still has to
      wait for the master process to let it complete.
      If a reader process does not have any read requests,
      it notifies the master that it completed matching and returns.*/
    DTF_DBG(VERBOSE_DBG_LEVEL, "Total %d rreqs and %d wreqs", fbuf->rreq_cnt, fbuf->wreq_cnt);
    if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
        send_ioreqs_by_var(fbuf, intracomp_io_flag);
    else if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
        send_ioreqs_by_block(fbuf, intracomp_io_flag);
    else 
		send_ioreqs_by_mst(fbuf, intracomp_io_flag);

//    t_part1 = MPI_Wtime() - t_part1;
//    t_part2 = MPI_Wtime();

	//TODO this flag is not really needed?
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
                if( (fbuf->writer_id == gl_my_comp_id) && fbuf->mst_info->is_master_flag  ){
                    do_matching(fbuf, intracomp_io_flag);
                    gl_stats.accum_do_matching_time += MPI_Wtime() - t_st;
                } else
                    gl_stats.idle_time += MPI_Wtime() - t_st;
            //}
    }

   
    reset_fbuf_match(fbuf);
    gl_stats.accum_match_time += MPI_Wtime() - t_start;

	
	{
		int is_scale = 0;
		
		char *c = getenv("DTF_SCALE");
		if(c != NULL)
			is_scale = atoi(c);
			
	    if(!is_scale){
			if(fbuf->mst_info->is_master_flag)
				clean_iodb(fbuf->mst_info->iodb, fbuf->nvars);
			
			delete_ioreqs(fbuf,finalize);
		}
	}

//    t_part2 = MPI_Wtime() - t_part2;
//    t_part3 = MPI_Wtime();
	
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);

    if(!intracomp_io_flag && (gl_my_comp_id == fbuf->writer_id)){
		MPI_Barrier(fbuf->comm);
		if(gl_my_rank == fbuf->root_writer){
			/*Notify reader writer finished this matching. Needed for multi-iterative cases*/
			dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_CONFIRM_TAG);
			int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, READ_DONE_CONFIRM_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
			CHECK_MPI(err);
			ENQUEUE_ITEM(msg, gl_msg_q);
			progress_msg_queue();
		}
	}
								
	DTF_DBG(VERBOSE_ERROR_LEVEL, "Stat: Time for matching %.4f", MPI_Wtime() - t_start);
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

        var = fbuf->vars[var_id];
		MPI_Type_size(var->dtype, &def_el_sz);
		/*Make sure that we can fit in at least one
		  element*/
		if( (var->ndims*sizeof(MPI_Offset)+1)*2 + def_el_sz+/*padding*/def_el_sz%sizeof(MPI_Offset) > sbufsz){
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: data message upper limit is \
					too small (current value %d). Please increase by setting up DFT_DATA_MSG_SIZE_LIMIT", gl_conf.data_msg_size_limit);
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
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
                int bt_bench = 0;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Will copy subblock (strt->cnt):");
                for(i=0; i < var->ndims; i++)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld \t orig: \t %lld\t --> %lld)", new_start[i], new_count[i], start[i], count[i]);

                /*Copy the block to send buffer*/

                char *s = getenv("BT_BENCH");
                if(s != NULL)  //HACK: don't check  for type mismatch for BT benchmark
                    bt_bench = 1;
                if(!bt_bench && (var->dtype != ioreq->dtype)){
                //if(var->dtype != ioreq->dtype){
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
                        if(bt_bench)
                            start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, new_start) * def_el_sz;
                        else
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
        var = fbuf->vars[var_id];
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
        //TODO move reqs to vars
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
        int bt_bench = 0;
         char *s = getenv("BT_BENCH");
        if(s != NULL)  //HACK: don't check  for type mismatch for BT benchmark
            bt_bench = 1;
        if(!bt_bench && (ioreq->dtype != var->dtype)){
        //if(ioreq->dtype != var->dtype){
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
        assert(nelems <= ioreq->user_buf_sz - ioreq->get_sz);

       DTF_DBG(VERBOSE_DBG_LEVEL, "Will get %d elems for var %d", nelems, var->id);
        if(var->ndims == 0){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Copy scalar variable");
            get_put_data(var, ioreq->dtype, ioreq->user_buf, NULL,
                            NULL, NULL, NULL, data, DTF_WRITE, type_mismatch);
        } else if(cnt_mismatch <= 1){
            MPI_Offset start_cpy_offt;
                if(bt_bench)
                    start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, start) * def_el_sz;
                else
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
            if(bt_bench)
                ioreq->get_sz += (MPI_Offset)(nelems*def_el_sz);
            else
            ioreq->get_sz += (MPI_Offset)(nelems*req_el_sz);


        DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, Got %d (expect %d)", ioreq->id, ioreq->var_id, (int)ioreq->get_sz, (int)ioreq->user_buf_sz);
        assert(ioreq->get_sz<=ioreq->user_buf_sz);
        if(ioreq->get_sz == ioreq->user_buf_sz){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Complete req %d (left %d), var %d, ", ioreq->id, fbuf->rreq_cnt-1, var->id);
            //for(i = 0; i < var->ndims; i++)
              //  DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", ioreq->start[i], ioreq->count[i]);
            //~ int nnels = 1;
			
			//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "------------COMPLETE IOREQ--------:");
			//~ for(i = 0; i < var->ndims; i++){
				//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "  %lld --> %lld", ioreq->start[i], ioreq->count[i]);
				//~ nnels *= ioreq->count[i];
			//~ }
			//~ for(i = 0; i < nnels; i++)
				//~ printf("%.2f\t", ((double*)ioreq->user_buf)[i]);
			//~ printf("\n");
			 
            if(gl_conf.do_checksum){
                double chsum = compute_checksum(ioreq->user_buf, var->ndims, ioreq->count, ioreq->dtype);
                DTF_DBG(VERBOSE_DBG_LEVEL, "Chsum for req is %.4f", chsum);
            }
            //delete this ioreq
            delete_ioreq(fbuf, &ioreq);

            if(fbuf->rreq_cnt == 0){
				int err;
				
                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE:Completed all rreqs for file %s", fbuf->file_path);
				
				dtf_msg_t *msg = new_dtf_msg(NULL, 0, READ_DONE_TAG);
				
				
				if(gl_my_comp_id == fbuf->writer_id)
					err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer,
                          READ_DONE_TAG, gl_comps[fbuf->writer_id].comm, &(msg->req));
                else
					err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader,
                          READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
                          
				CHECK_MPI(err);
				
				//Reader ranks except for root complete the send immediately
				if(gl_my_comp_id == fbuf->reader_id && gl_my_rank != fbuf->root_reader){
					err = MPI_Wait(&(msg->req), MPI_STATUS_IGNORE);
					CHECK_MPI(err);
					dtf_free(msg, sizeof(dtf_msg_t));
					fbuf->done_matching_flag = 1;
				} else 			
					ENQUEUE_ITEM(msg, gl_msg_q);
            } 
            //else {
				//io_req_t *tmp = fbuf->ioreqs;
				//DTF_DBG(VERBOSE_DBG_LEVEL, "Left rreqs:");
				//while(tmp != NULL){
                    //DTF_DBG(VERBOSE_DBG_LEVEL, "%u (varid %d)", tmp->id, tmp->var_id);
                    //tmp = tmp->next;
				//}
			//}
        }
    } //while(bufsz)

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
        t_st = MPI_Wtime();
        err = MPI_Test(&(msg->req), &flag, &status);
        CHECK_MPI(err);
        if(flag){
            gl_stats.accum_comm_time += MPI_Wtime() - t_st;
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
        } else{
            gl_stats.idle_time += MPI_Wtime() - t_st;
            msg = msg->next;
        }
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
        //~ case IO_CLOSE_FILE_TAG:
            //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_CLOSE_FILE_TAG");
            //~ break;
        //~ case IO_OPEN_FILE_FLAG:
            //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_OPEN_FILE_FLAG");
            //~ break;
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
        case READ_DONE_CONFIRM_TAG:
			DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag READ_DONE_CONFIRM_TAG");
			break;
        case SKIP_MATCH_TAG:
			DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag SKIP_MATCH_TAG");
			break;
        default:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag unknown %d", tag);
            assert(0);
    }
}

/*Check if there are any unfinished sends of ioreqs and cancel them*/
void cancel_send_ioreqs()
{
	dtf_msg_t *msg, *tmp;
    int flag, err;
    MPI_Status status;

    if(gl_msg_q == NULL)
       return;
    msg = gl_msg_q;
    while(msg != NULL){
        if(msg->tag != IO_REQS_TAG){
			msg = msg->next;
			continue;
		}
		DTF_DBG(VERBOSE_ERROR_LEVEL, "Try to cancel msg %p", (void*)msg);
		//try to cancel the request
		err = MPI_Cancel(&(msg->req));
		CHECK_MPI(err);
		err = MPI_Wait(&(msg->req), &status); 
		CHECK_MPI(err);
	    err = MPI_Test_cancelled( &status, &flag );
		if (!flag) 
			DTF_DBG(VERBOSE_ERROR_LEVEL," DTF Warning: Failed to cancel an Isend request\n" );
		else 
			DTF_DBG(VERBOSE_ERROR_LEVEL, "Canceled");
			
		tmp = msg->next;
		DEQUEUE_ITEM(msg, gl_msg_q);		
		delete_dtf_msg(msg);
		msg = tmp;
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
            t_st = MPI_Wtime();
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                gl_stats.idle_time += MPI_Wtime() - t_st;
                break;
            }
            src = status.MPI_SOURCE;

            print_recv_msg(status.MPI_TAG);

            switch(status.MPI_TAG){
                case FILE_INFO_REQ_TAG:
					//TODO remove this tag and use OPEN_FILE for everything
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv file info req from rank %d (comp %s)", src, gl_comps[comp].name);
                
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_INFO_REQ_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
					fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
					assert(fbuf != NULL);
                    assert(fbuf->root_writer == gl_my_rank);
                    
                    
					DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer, process the file info request");
					fbuf->root_reader = src;
					if(fbuf->iomode == DTF_IO_MODE_FILE){
						fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;
						if(fbuf->is_ready) //writer has already closed the file
							notify_file_ready(fbuf);
					} else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
						send_file_info(fbuf, fbuf->root_reader);
					}					
                    break;
                case IO_REQS_TAG:
                    t_st = MPI_Wtime();
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
                    gl_stats.master_time += MPI_Wtime() - t_st;
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
                case SKIP_MATCH_TAG:
					t_st = MPI_Wtime();
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, SKIP_MATCH_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    
                    if(gl_my_rank == fbuf->root_writer){
						DTF_DBG(VERBOSE_DBG_LEVEL, "Notify masters & workgroup that matching is finished (reader skipped)");
						notify_masters(fbuf, SKIP_MATCH_TAG);
						notify_workgroup(fbuf, MATCH_DONE_TAG);
						fbuf->done_matching_flag = 1;
					} else {
						assert(fbuf->mst_info->is_master_flag);
						notify_workgroup(fbuf, MATCH_DONE_TAG);
						fbuf->done_matching_flag = 1;
					}
                    gl_stats.master_time += MPI_Wtime() - t_st;
					break;
                case READ_DONE_TAG:
                    t_st = MPI_Wtime();
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, READ_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
					
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv read done for file %s from %d in comp %d (tot %d)", fbuf->file_path, 
							src, comp, fbuf->mst_info->nrranks_completed);
					assert(fbuf->mst_info->nrranks_completed < fbuf->mst_info->nranks_opened);
					
					fbuf->mst_info->nrranks_completed++;
										
					if( (comp == fbuf->writer_id) && (fbuf->mst_info->nrranks_completed == fbuf->mst_info->nranks_opened)){ 
						//finished intracomp matching
						assert(gl_my_rank == fbuf->root_writer);
						DTF_DBG(VERBOSE_DBG_LEVEL, "Notify masters & workgroup that intracomp matching is finished");
						notify_masters(fbuf, MATCH_DONE_TAG);
						notify_workgroup(fbuf, MATCH_DONE_TAG);
						fbuf->done_matching_flag = 1;
						DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE, Done matching at %.3f", MPI_Wtime()-gl_stats.walltime);		
					} else if(comp == fbuf->reader_id){
						
						if(gl_my_comp_id == fbuf->writer_id){
							
							assert(fbuf->mst_info->is_master_flag);
							notify_workgroup(fbuf, MATCH_DONE_TAG);
							fbuf->done_matching_flag = 1;
							DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
							
						} else if((fbuf->root_reader == gl_my_rank) && (fbuf->mst_info->nrranks_completed == fbuf->mst_info->nranks_opened)){
								int i;
								DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer masters readers completed matching");
								for(i = 0; i < fbuf->mst_info->nmasters; i++){
									int err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->mst_info->masters[i], READ_DONE_TAG, gl_comps[fbuf->writer_id].comm);
									CHECK_MPI(err);
								}
								fbuf->done_matching_flag = 1;
								DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
						}
						
					}
				
                    gl_stats.master_time += MPI_Wtime() - t_st;
                    break;
                case MATCH_DONE_TAG:

                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, MATCH_DONE_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    //if(fbuf->mst_info->is_master_flag){
                        //t_st = MPI_Wtime();
                        //DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that req matching for %s completed", fbuf->file_path);
                        ///*Tell other writer ranks that they can complete matching*/
                        //notify_workgroup(fbuf, MATCH_DONE_TAG);
                        //gl_stats.master_time += MPI_Wtime() - t_st;
                    //}
                    
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
                    fbuf->done_matching_flag = 1;
                    if(gl_conf.iodb_build_mode == IODB_BUILD_RANK)
						cancel_send_ioreqs();
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
                case READ_DONE_CONFIRM_TAG:
					err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, READ_DONE_CONFIRM_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    assert(fbuf != NULL);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Recv done match confirm tag for %s from %d", fbuf->file_path, src);
                    fbuf->done_match_confirm_flag = 1;
					break;
                //~ case IO_CLOSE_FILE_TAG:
                    //~ err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, IO_CLOSE_FILE_TAG, gl_comps[comp].comm, &status);
                    //~ CHECK_MPI(err);
                    //~ fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    //~ assert(fbuf != NULL);
                    //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Recv close file tag for %s from %d", fbuf->file_path, src);
                    //~ //assert(!fbuf->rdr_closed_flag); //check that haven't received that before

                    //~ if(gl_my_rank == fbuf->root_writer){
                        //~ //DTF_DBG(VERBOSE_DBG_LEVEL, "Notify other masters that readers are closing the file");
                        //~ //notify_masters(fbuf, IO_CLOSE_FILE_TAG);

                        //~ //TODO for multi-cycle version will need to figure out with these states and notifications
                        //~ //note: this should be set when opening or creating the file!
                        //~ //if(fbuf->iomode == DTF_IO_MODE_FILE)
                          //~ //  fbuf->fready_notify_flag = RDR_HASNT_OPENED;
                    //~ }
                    //~ if(fbuf->mst_info->is_master_flag){
                        //~ t_st = MPI_Wtime();
                        //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that they can close the file %s", fbuf->file_path);
                        //~ notify_workgroup(fbuf, IO_CLOSE_FILE_TAG);
                        //~ gl_stats.master_time += MPI_Wtime() - t_st;
                    //~ }
                
                    //~ fbuf->rdr_closed_flag = 1;

                    //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Close flag set for file %s", fbuf->file_path);
//~ //                    if(fbuf->done_matching_flag){
//~ //                        if(fbuf->mst_info->is_master_flag){
//~ //                            /*Check that I don't have any read reqs incompleted*/
//~ //                            assert(fbuf->mst_info->iodb->nritems == 0);
//~ //                            /*Clean my iodb*/
//~ //                            clean_iodb(fbuf->mst_info->iodb);
//~ //                        }
//~ //                         /*Can delete write requests only if all ps-s have finished
//~ //                         matching.*/
//~ //                        delete_ioreqs(fbuf);
//~ //                    }
                    //~ break;
                case ROOT_MST_TAG:
                    src = status.MPI_SOURCE;
                    err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, ROOT_MST_TAG, gl_comps[comp].comm, &status);
                    CHECK_MPI(err);
                    DTF_DBG(VERBOSE_DBG_LEVEL,   "Receive ROOT_MST_TAG notif for %s", filename);
					
					{
						file_info_t *finfo = dtf_malloc(sizeof(file_info_t));
						assert(finfo != NULL);
						strcpy(finfo->filename, filename);
						finfo->root_writer = src;
						finfo->prev = NULL;
						finfo->next = NULL;
						if(gl_finfo_list == NULL)
							gl_finfo_list = finfo;
						else{ 
							finfo->next = gl_finfo_list;
							finfo->next->prev = finfo;
							gl_finfo_list = finfo;
						}
					}
					
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
                        t_st = MPI_Wtime();
                        DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that multiple matching for %s completed", fbuf->file_path);
                        /*Tell other writer ranks that they can complete matching*/
                        notify_workgroup(fbuf, DONE_MULTIPLE_FLAG);
                        gl_stats.master_time += MPI_Wtime() - t_st;
                    }
                    fbuf->done_matching_flag = 1;
                    fbuf->done_match_multiple_flag = 1;
                    break;
                //~ case IO_OPEN_FILE_FLAG:
                    //~ err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, IO_OPEN_FILE_FLAG, gl_comps[comp].comm, &status);
                    //~ CHECK_MPI(err);
                    //~ fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
                    //~ assert(fbuf != NULL);
                    //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Recv open file tag for %s from %d", fbuf->file_path, src);
                    //~ fbuf->rdr_closed_flag = 0;
                    //~ if(fbuf->mst_info->is_master_flag){
                        //~ t_st = MPI_Wtime();
                        //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Notify workgroup");
                        //~ notify_workgroup(fbuf, IO_OPEN_FILE_FLAG);

                        //~ if(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->root_writer == gl_my_rank)
							//~ fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;									
						
                        //~ gl_stats.master_time += MPI_Wtime() - t_st;
                    //~ }
                    //~ break;
                default:
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: unknown tag %d", status.MPI_TAG);
                    assert(0);
            }
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

    mst_info->nrranks_completed = 0;
    mst_info->iodb = NULL;
    return 0;
}

