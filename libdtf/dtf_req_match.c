#include "dtf_req_match.h"
#include "dtf_util.h"
#include "dtf.h"
#include <unistd.h>
#include "dtf_io_pattern.h"

static void shift_coord(int ndims, const MPI_Offset *bl_start,
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
    DTF_DBG(VERBOSE_ALL_LEVEL, "New start before adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_ALL_LEVEL, "\t %lld", subbl_start[i]);

    for(i = ndims - 1; i > 0; i--)
        if(subbl_start[i] == bl_start[i] + bl_count[i]){
            subbl_start[i] = bl_start[i];
            if( (subbl_start[i-1] != bl_start[i-1] + bl_count[i-1]) && (subbl_count[i-1] == 1)){
                subbl_start[i-1]++;
            }
        } else
            break;

    DTF_DBG(VERBOSE_ALL_LEVEL, "New start after adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_ALL_LEVEL, "\t %lld", subbl_start[i]);


//    DTF_DBG(VERBOSE_DBG_LEVEL, "Copied subblock. Shift start:");
//    for(i = 0; i < var->ndims; i++)
//        DTF_DBG(VERBOSE_DBG_LEVEL, "   %lld\t -->\t %lld", bl_start[i], subbl_start[i]);
}

void delete_ioreq(file_buffer_t *fbuf, int varid, io_req_t **ioreq)
{
	DTF_DBG(VERBOSE_ALL_LEVEL, "Delete req %d, cur wreqs %d, rreqs %d", (*ioreq)->id,
			fbuf->wreq_cnt, fbuf->rreq_cnt);

    dtf_var_t *var = fbuf->vars[varid];

    if( (*ioreq)->is_buffered)
        dtf_free((*ioreq)->user_buf, (size_t) (*ioreq)->user_buf_sz);

    if(*ioreq == var->ioreqs)
        var->ioreqs = var->ioreqs->next;
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
    int varid;
    dtf_var_t *var;
    DTF_DBG(VERBOSE_ALL_LEVEL, "Delete io requests for file %s", fbuf->file_path);
	for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			if(ioreq->is_permanent && !finalize){
				ioreq = ioreq->next;
				continue;
			}
			tmp = ioreq->next;
			delete_ioreq(fbuf, varid, &ioreq);
			ioreq = tmp;
		}
		if(finalize)
			assert(var->ioreqs == NULL);
	}

    if(finalize){
		assert(fbuf->rreq_cnt == 0);
		assert(fbuf->wreq_cnt == 0);
	}
}




/*match subblocks of data*/
static void do_matching(file_buffer_t *fbuf)
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
    if(!fbuf->my_mst_info->iodb->updated_flag){ //no new info since last time matching was done, ignore
        return;
	}
    if(fbuf->my_mst_info->iodb->witems == NULL || fbuf->my_mst_info->iodb->ritems == NULL){
		return;
	}
    fbuf->my_mst_info->iodb->updated_flag = 0; //reset

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
    ritem = fbuf->my_mst_info->iodb->ritems;
    while(ritem != NULL){

        t_st = MPI_Wtime();

        n_matched_blocks = 0;

        for(i = 0; i < allocd_nwriters; i++){
            sbuf[i] = NULL;
            bufsz[i] = 0;
            offt[i] = 0;
            //writers[i] = -1;
        }

        nwriters = 0;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Matching rreq from rank %d", ritem->rank);
        //~ if(ritem->comm == gl_comps[gl_my_comp_id].comm)
            //~ DTF_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in my comp", ritem->rank);
        //~ else
            //~ DTF_DBG(VERBOSE_ALL_LEVEL, "rreq from rank %d in other comp", ritem->rank);


        rblock = ritem->dblocks;
        while(rblock != NULL){
            ntimes_while++;
            var_id = rblock->var_id;
            t_st = MPI_Wtime();
			assert(fbuf->my_mst_info->iodb->witems != NULL);
			
            witem = fbuf->my_mst_info->iodb->witems[var_id];

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
                    DTF_DBG(VERBOSE_DBG_LEVEL, "didnt find block for var %d", var_id);
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
                      and reader rank*/
                    if(offt[rank_idx] == 0) {
                    if(bufsz[rank_idx] < MAX_FILE_NAME + sizeof(MPI_Offset)){
                            //extend
                            unsigned char *tmp;
                            tmp = realloc(sbuf[rank_idx], bufsz[rank_idx] + MAX_FILE_NAME + mlc_buf);
                            assert(tmp != NULL);
                            sbuf[rank_idx] = tmp;
                            bufsz[rank_idx] += MAX_FILE_NAME + mlc_buf;
                        }

                        memcpy(sbuf[rank_idx], fbuf->file_path, MAX_FILE_NAME);
                        offt[rank_idx] += MAX_FILE_NAME ;
						*(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)ritem->rank; //to whom writer should send the data
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
                    dtf_msg_t *msg = new_dtf_msg(sbuf[i], bufsz[i], DTF_UNDEFINED, IO_DATA_REQ_TAG);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Send data req to wrt %d", writers[i]);
                    err = MPI_Isend((void*)sbuf[i], offt[i], MPI_BYTE, writers[i], IO_DATA_REQ_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
                    CHECK_MPI(err);
                    ENQUEUE_ITEM(msg, gl_out_msg_q);
                }
            }

            if(my_idx != -1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Parse data req to myself");
                send_data(fbuf, sbuf[my_idx] + MAX_FILE_NAME, offt[my_idx] - MAX_FILE_NAME);
                dtf_free(sbuf[my_idx], bufsz[my_idx]);
            }
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for rank %d", ritem->rank);
        }

        /*If we matched all chunks for this rank, then delete this ritem*/
        if(ritem->nblocks == 0){
            read_db_item_t *tmp = ritem;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched all. Delete ritem of rank %d (left ritems %d). ", ritem->rank,  (int)(fbuf->my_mst_info->iodb->nritems - 1));
            if(ritem == fbuf->my_mst_info->iodb->ritems)
                fbuf->my_mst_info->iodb->ritems = ritem->next;
            if(ritem->prev != NULL)
                ritem->prev->next = ritem->next;
            if(ritem->next != NULL)
                ritem->next->prev = ritem->prev;
            ritem = ritem->next;
            dtf_free(tmp, sizeof(read_db_item_t));
            fbuf->my_mst_info->iodb->nritems--;
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
	gl_stats.t_do_match += MPI_Wtime() - t_start;

    DTF_DBG(VERBOSE_DBG_LEVEL, "after matching: %d ritems", (int)fbuf->my_mst_info->iodb->nritems);
}


static void parse_ioreqs(file_buffer_t *fbuf, void *buf, int bufsz, int rank, MPI_Comm comm)
{
    int var_id, rw_flag, i;
    dtf_var_t *var = NULL;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
	read_db_item_t *ritem = NULL;

    double t_st = MPI_Wtime();
    
   // offt += MAX_FILE_NAME;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Start parsing reqs for file %s", fbuf->file_path);
    if(comm == gl_comps[fbuf->reader_id].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Reqs are from reader");
    else {
        assert(comm == gl_comps[fbuf->writer_id].comm);
        DTF_DBG(VERBOSE_DBG_LEVEL, "Req are from writer");
    }

	if(fbuf->done_matching_flag){
		DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: received io reqs for file %s after matching already finished, discard ioreqs", fbuf->file_path);
		return;
	}

    while(offt != (size_t)bufsz){
        rw_flag = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        gl_stats.iodb_nioreqs++;
        assert(var_id < fbuf->nvars);
		var = fbuf->vars[var_id];

        if(rw_flag == DTF_READ){

			 if(ritem == NULL){

				ritem = fbuf->my_mst_info->iodb->ritems;
				while(ritem != NULL){
					if(ritem->rank == rank)
						break;
					ritem = ritem->next;
				}

				if(ritem == NULL){
					ritem = (read_db_item_t*)dtf_malloc(sizeof(read_db_item_t));
					assert(ritem != NULL);
					ritem->rank = rank;
					ritem->dblocks = NULL;
					ritem->nblocks = 0;
					ritem->last_block = NULL;
					ritem->next = NULL;
					ritem->prev = NULL;
					//enqueue
					if(fbuf->my_mst_info->iodb->ritems == NULL)
						fbuf->my_mst_info->iodb->ritems = ritem;
					else{
						ritem->next = fbuf->my_mst_info->iodb->ritems;
						fbuf->my_mst_info->iodb->ritems->prev = ritem;
						fbuf->my_mst_info->iodb->ritems = ritem;
					}
					fbuf->my_mst_info->iodb->nritems++;
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
            
            DTF_DBG(VERBOSE_DBG_LEVEL, "Added dblock:");
			for(i = 0; i < var->ndims; i++)
				DTF_DBG(VERBOSE_DBG_LEVEL, "%lld  ->  %lld", dblock->start[i], dblock->count[i]);
			
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
            DTF_DBG(VERBOSE_DBG_LEVEL, "ritems %d, cur item (r %d) %lld blocks",(int)fbuf->my_mst_info->iodb->nritems, ritem->rank, ritem->nblocks);
        } else { /*DTF_WRITE*/
			int i;
            write_db_item_t *witem;
            /*Allow write requests only from the writer*/
            assert(comm == gl_comps[fbuf->writer_id].comm);
            
            if(fbuf->my_mst_info->iodb->witems == NULL){
				fbuf->my_mst_info->iodb->witems = dtf_malloc(fbuf->nvars*sizeof(write_db_item_t*));
				assert(fbuf->my_mst_info->iodb->witems != NULL);
				for(i = 0; i < fbuf->nvars; i++)
					fbuf->my_mst_info->iodb->witems[i] = NULL;
			}
            if(fbuf->my_mst_info->iodb->witems[var_id] == NULL){
                fbuf->my_mst_info->iodb->witems[var_id] = (write_db_item_t*)dtf_malloc(sizeof(write_db_item_t));
                assert(fbuf->my_mst_info->iodb->witems[var_id] != NULL);
                witem = fbuf->my_mst_info->iodb->witems[var_id];

                if(var->ndims > 0)
                    witem->dblocks = (void*)RBTreeCreateBlocks(rb_key_cmp, NullFunction, rb_destroy_node_info, rb_print_key, rb_print_info, 0);
                else{
                    witem->dblocks = dtf_malloc(sizeof(block_t));
                    assert(witem->dblocks != NULL);
				}
                witem->nblocks = 0;
                witem->ndims = var->ndims;
            } else
				witem = fbuf->my_mst_info->iodb->witems[var_id];
				
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
                //~ DTF_DBG(VERBOSE_ALL_LEVEL, "Insert block for var %d", var_id);
               
				//~ for(i = 0; i < var->ndims; i++)
					//~ DTF_DBG(VERBOSE_DBG_LEVEL, "%lld  ->  %lld", info->blck->start[i], info->blck->count[i]);
				
				rb_red_blk_node *bl_node = RBTreeInsertBlock(witem->dblocks, info);
				assert(bl_node != NULL);
				dtf_free(info, sizeof(insert_info));

			} else{
                ((block_t*)witem->dblocks)->start = NULL;
                ((block_t*)witem->dblocks)->count = NULL;
                ((block_t*)witem->dblocks)->rank = rank;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Inserted scalar var");
            }
            witem->nblocks++;
        }
    }
    assert(offt == (size_t)bufsz);
	fbuf->my_mst_info->iodb->updated_flag = 1;
    gl_stats.parse_ioreq_time += MPI_Wtime() - t_st;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished parsing reqs. (mem %lu)", gl_stats.malloc_size);
    
    do_matching(fbuf);
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
    } else {
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
    ioreq->get_sz = 0;
    ioreq->prev = NULL;
    ioreq->dtype = dtype;
    ioreq->rw_flag = rw_flag;
    ioreq->is_permanent = 0;

    if( (rw_flag == DTF_WRITE) && gl_conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
        ioreq->checksum = compute_checksum(buf, ndims, count, dtype);
        DTF_DBG(VERBOSE_DBG_LEVEL, "chsum for req %d is %.4f", ioreq->id, ioreq->checksum);
    } else
        ioreq->checksum = 0;
    gl_stats.nioreqs++;
    return ioreq;
}

/*The metadata is distributed among masters based on the var id.*/
void send_ioreqs_by_var(file_buffer_t *fbuf)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int mst = 0, varid;
    int nrreqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    int data_to_send = 0;
    int nmasters;
    int idx, err;
    master_info_t *mst_info;

    if(fbuf->rreq_cnt == 0 && fbuf->wreq_cnt == 0)
        return;

    if(gl_my_comp_id == fbuf->writer_id)
        mst_info = fbuf->my_mst_info;
    else
        mst_info = fbuf->cpl_mst_info;
    nmasters = mst_info->nmasters;

    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    assert(sbuf != NULL);
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(offt != NULL);
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(bufsz != NULL);

    /*Distribute ioreqs between the masters based on var id*/

    //alloc mem
    for(mst = 0; mst < nmasters; mst++){
        bufsz[mst] = 0;
        sbuf[mst] = NULL;
    }

    nrreqs = 0;
    for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			
			if(ioreq->rw_flag == DTF_READ)
				nrreqs++;
				
			if(ioreq->sent_flag){
				//All the following reqs should have been sent already
				break;
			}
			data_to_send = 1;

			mst = varid % nmasters;
			bufsz[mst] += sizeof(MPI_Offset)*2 + var->ndims*2*sizeof(MPI_Offset);
			ioreq = ioreq->next;
		}
	}

    if(fbuf->reader_id == gl_my_comp_id && nrreqs == 0){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
		int flag; 
		
		dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, READ_DONE_TAG);
		err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
				  READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
		CHECK_MPI(err);
		err = MPI_Test(&(msg->req), &flag, MPI_STATUS_IGNORE);
		CHECK_MPI(err);
		ENQUEUE_ITEM(msg, gl_out_msg_q);
		
		if(gl_my_rank != fbuf->root_reader){
			fbuf->done_matching_flag = 1;
		}	
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    for(mst = 0; mst < nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;

        bufsz[mst] += MAX_FILE_NAME ;
        DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu for mst %d", bufsz[mst], mst);
        sbuf[mst] = dtf_malloc(bufsz[mst]);
        assert(sbuf[mst] != NULL);

        memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
        offt[mst] = MAX_FILE_NAME ;
    }

    for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			if(ioreq->sent_flag){
				//All the following reqs should have been sent already
				break;
			}
			mst = varid % nmasters;
			DTF_DBG(VERBOSE_DBG_LEVEL, "Will send ioreq %d (varid %d) to mst %d", ioreq->id, varid, mst); 
			/*Store var_id, rw_flag, start[] and count[]*/
			*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
			offt[mst] += sizeof(MPI_Offset);
			*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)varid;
			offt[mst] += sizeof(MPI_Offset);
			memcpy(sbuf[mst]+offt[mst], ioreq->start, var->ndims*sizeof(MPI_Offset));
			offt[mst] += var->ndims*sizeof(MPI_Offset);
			memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));
			offt[mst] += var->ndims*sizeof(MPI_Offset);
			ioreq->sent_flag = 1;
			ioreq = ioreq->next;
		}
	}

    idx = -1;
    for(mst = 0; mst < nmasters; mst++){
        if(bufsz[mst] == 0){
            continue;
        }
        if( (fbuf->writer_id == gl_my_comp_id) && (mst_info->masters[mst] == gl_my_rank)){
            idx = mst;
        } else {
            int flag;
            MPI_Status status;
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], DTF_UNDEFINED, IO_REQS_TAG);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu), comp id %d (mine %d)", mst_info->masters[mst], offt[mst], bufsz[mst], fbuf->writer_id, gl_my_comp_id);
            assert(offt[mst] == bufsz[mst]);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, mst_info->masters[mst], IO_REQS_TAG,
                            gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            err = MPI_Test(&(msg->req), &flag, &status);
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_out_msg_q);
        }
    }
    
    if(idx != -1){
		assert(offt[idx] == bufsz[idx]);
        parse_ioreqs(fbuf,sbuf[idx], offt[idx], gl_my_rank, gl_comps[gl_my_comp_id].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    } 
    
    fbuf->t_last_sent_ioreqs = MPI_Wtime();
	fbuf->has_unsent_ioreqs = 0;
fn_exit:
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz,nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}


/*This version divides the variable data among masters along the var's biggest dimension.
  If there is an unlimited dimension, the division happens along it*/
void send_ioreqs_by_block(file_buffer_t *fbuf)
{
    dtf_var_t *var = NULL;
    io_req_t *ioreq;
    int varid;
    int mst = 0;    //only one master for now
    int nrreqs = 0;
    unsigned char **sbuf;
    size_t *bufsz, *offt;
    size_t mlc_chunk = 512*1024;
    int data_to_send = 0;
    int nmasters;
    int idx, err, i;
    MPI_Offset block_range, shift;
    MPI_Offset *strt, *cnt;
    unsigned block_cnt = 0;
    master_info_t *mst_info;

    if(fbuf->rreq_cnt == 0 && fbuf->wreq_cnt==0)
        return;

    if(gl_my_comp_id == fbuf->writer_id)
        mst_info = fbuf->my_mst_info;
    else
        mst_info = fbuf->cpl_mst_info;

    nmasters = mst_info->nmasters;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Sending I/O reqs by block");

    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    assert(sbuf != NULL);
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(offt != NULL);
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    assert(bufsz != NULL);

    /*Distribute ioreqs between the masters based on the first dimension of the var*/

    for(mst = 0; mst < nmasters; mst++){
        bufsz[mst] = 0;
        offt[mst] = 0;
        sbuf[mst] = NULL;
    }

    nrreqs = 0;
	for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			if(ioreq->rw_flag == DTF_READ)
				nrreqs++;
				
			if(ioreq->sent_flag){
				//All the following reqs should have been sent already
				break;
			}
			data_to_send = 1;
			
			ioreq = ioreq->next;
		}
	}
   if(fbuf->reader_id == gl_my_comp_id && nrreqs == 0){
		 DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
         int flag = 0;
		
		dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, READ_DONE_TAG);
		err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
				  READ_DONE_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
		CHECK_MPI(err);
		err = MPI_Test(&(msg->req), &flag, MPI_STATUS_IGNORE);
		CHECK_MPI(err);
		ENQUEUE_ITEM(msg, gl_out_msg_q);
		
		if(gl_my_rank != fbuf->root_reader){
			fbuf->done_matching_flag = 1;
		}
		
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			if(ioreq->sent_flag){
				break;
			}
			DTF_DBG(VERBOSE_ALL_LEVEL, "Will send ioreq for var %d:", varid);
			for(i = 0; i < var->ndims; i++)
				DTF_DBG(VERBOSE_ALL_LEVEL, "%lld  ->  %lld", ioreq->start[i], ioreq->count[i]);

			if(var->ndims == 0){
				mst = varid % nmasters;
				if(bufsz[mst] == 0){
						sbuf[mst] = dtf_malloc(MAX_FILE_NAME + mlc_chunk);
						assert(sbuf[mst] != NULL);
						bufsz[mst] = MAX_FILE_NAME + mlc_chunk;

						memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
						offt[mst] = MAX_FILE_NAME ;
				}
				/*Only save the rw flag and the var id*/
				*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
				offt[mst] += sizeof(MPI_Offset);
				*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)varid;
				offt[mst] += sizeof(MPI_Offset);
				block_cnt++;
			} else {

				if(var->shape[0] == DTF_UNLIMITED)
					block_range = DEFAULT_BLOCK_SZ_RANGE;
				else if(gl_conf.iodb_range > 0)
					block_range = gl_conf.iodb_range;
				else{
					block_range = (MPI_Offset)(var->shape[0]/nmasters);
					if(var->shape[0]%nmasters>0)
						block_range++;
				}
				if(block_range == 0)
					block_range = DEFAULT_BLOCK_SZ_RANGE;

				DTF_DBG(VERBOSE_ALL_LEVEL, "Var %d, block_range %lld", var->id, block_range);

				shift = 0;
				while(ioreq->start[0] + shift < ioreq->start[0] + ioreq->count[0]){
					block_cnt++;
					if( (block_range == DEFAULT_BLOCK_SZ_RANGE) || (block_range == gl_conf.iodb_range) )
						mst = (int)( ((ioreq->start[0] + shift)/block_range) % nmasters);
					else
						mst = (int)((ioreq->start[0] + shift)/block_range);

					DTF_DBG(VERBOSE_ALL_LEVEL, "mst %d (%lld)", mst, ioreq->start[0] + shift);
					assert(mst < nmasters);
					if(bufsz[mst] == 0){
						sbuf[mst] = dtf_malloc(MAX_FILE_NAME + mlc_chunk);
						assert(sbuf[mst] != NULL);
						bufsz[mst] = MAX_FILE_NAME + mlc_chunk;

						memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
						offt[mst] = MAX_FILE_NAME ;
					}

					if(var->ndims*sizeof(MPI_Offset)*3+offt[mst] > bufsz[mst]){
						size_t ext_sz = mlc_chunk;
						while(bufsz[mst]+ext_sz < var->ndims*sizeof(MPI_Offset)*3+offt[mst])
							ext_sz += mlc_chunk;
						DTF_DBG(VERBOSE_ALL_LEVEL, "bufsz %lu, ext sz %lu", bufsz[mst], ext_sz );
						
						void *tmp = realloc(sbuf[mst], bufsz[mst] + ext_sz);
						assert(tmp != NULL);
						gl_stats.malloc_size += ext_sz;
						bufsz[mst] += ext_sz;
					}

					/*Store var_id, rw_flag, start[] and count[]*/
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
					offt[mst] += sizeof(MPI_Offset);
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)varid;
					offt[mst] += sizeof(MPI_Offset);
					memcpy(sbuf[mst]+offt[mst], ioreq->start, var->ndims*sizeof(MPI_Offset));
					/*Adjust corresponding start coordinate*/
					*(MPI_Offset*)(sbuf[mst]+offt[mst] + 0*sizeof(MPI_Offset)) = ioreq->start[0] + shift;
					strt = (MPI_Offset*)(sbuf[mst]+offt[mst]);
					offt[mst] += var->ndims*sizeof(MPI_Offset);
					
					memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));
					cnt = (MPI_Offset*)(sbuf[mst]+offt[mst]);
					/*Adjust corresponding count*/
					//~ if(ioreq->count[0] - shift > block_range - (ioreq->start[0]+shift-block_range*mst) ){
						//~ *(MPI_Offset*)(sbuf[mst]+offt[mst] + 0*sizeof(MPI_Offset)) =  block_range - (ioreq->start[0]+shift-block_range*mst);
						//~ shift +=  block_range - (ioreq->start[0]+shift-block_range*mst);
					//~ } else {
						//~ cnt[0] = ioreq->count[0] - shift;
						//~ shift = ioreq->count[0];
					//~ }

					if(ioreq->start[0]+cnt[0] > ioreq->start[0] + shift - (ioreq->start[0] + shift)%block_range + block_range) {
						cnt[0] = block_range -  (ioreq->start[0] + shift)%block_range;

					} else {
						cnt[0] = ioreq->count[0] - shift;

					}
					shift += cnt[0];
					offt[mst] += var->ndims*sizeof(MPI_Offset);
					
					DTF_DBG(VERBOSE_ALL_LEVEL, "Will send info to mst %d about block (shift along %d by %lld):", mst, 0, shift);
					for(i = 0; i < var->ndims; i++)
						DTF_DBG(VERBOSE_ALL_LEVEL, "%lld  ->  %lld", strt[i], cnt[i]);

				}
				ioreq->sent_flag = 1;
				DTF_DBG(VERBOSE_ALL_LEVEL, "Packed %u blocks", block_cnt);
			}
			ioreq = ioreq->next;
		}
	}
    idx = -1;

    for(mst = 0; mst < nmasters; mst++){
        if(bufsz[mst] == 0){
            continue;
        }
        if( (fbuf->writer_id == gl_my_comp_id) && (mst_info->masters[mst] == gl_my_rank)){
            idx = mst;
        } else {
            int flag;
            MPI_Status status;
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], DTF_UNDEFINED, IO_REQS_TAG);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu), comp id %d (mine %d)", mst_info->masters[mst], 		offt[mst], bufsz[mst], fbuf->writer_id, gl_my_comp_id);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, mst_info->masters[mst], IO_REQS_TAG,
                            gl_comps[fbuf->writer_id].comm, &(msg->req));
            CHECK_MPI(err);
            err = MPI_Test(&(msg->req), &flag, &status);
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_out_msg_q);
        }
    }

    //gl_stats.t_comm += MPI_Wtime() - t_start_comm;

    if(idx != -1){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Parse reqs to myself: %d", (int)offt[idx] - MAX_FILE_NAME);
        parse_ioreqs(fbuf, (unsigned char*)sbuf[idx] + MAX_FILE_NAME, offt[idx] - MAX_FILE_NAME, gl_my_rank, gl_comps[gl_my_comp_id].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    }
    
    fbuf->t_last_sent_ioreqs = MPI_Wtime();
    fbuf->has_unsent_ioreqs = 0;

fn_exit:
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz, nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}

void match_ioreqs_multiple()
{
	file_buffer_t *fbuf;
	int compl_file_cnt = 0;
	int file_cnt = 0;
	double t_start = MPI_Wtime();
	
	/*First check for how many files need to complete transfer 
	 * and send any unsent I/O requests*/
	fbuf = gl_filebuf_list;
	while(fbuf != NULL){
		if(fbuf->is_transfering){ 
			file_cnt++;
			DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is in active transfer", fbuf->file_path);
			
			if(strstr(fbuf->file_path, "hist.d")!=NULL)
				gl_stats.st_mtch_hist = MPI_Wtime()-gl_stats.walltime;
			else if(strstr(fbuf->file_path, "anal.d")!=NULL)
				gl_stats.st_mtch_rest = MPI_Wtime()-gl_stats.walltime;
			
			if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
				send_ioreqs_by_var(fbuf);
			else //if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
				send_ioreqs_by_block(fbuf);
		}
		fbuf = fbuf->next;
	}

	progress_comm();
	
	while(compl_file_cnt != file_cnt){
		
		fbuf = gl_filebuf_list;
		while(fbuf != NULL){
			if(fbuf->is_transfering && fbuf->done_matching_flag){
				
				compl_file_cnt++;
				
				/*Reset flags*/
				if (fbuf->writer_id == gl_my_comp_id)
					/*Reset flag for future matchings*/
					fbuf->my_mst_info->iodb->updated_flag = 0;
				if(fbuf->my_mst_info->my_mst == gl_my_rank)
					fbuf->my_mst_info->nread_completed = 0;
				fbuf->done_matching_flag = 0;
				fbuf->is_transfering = 0;
				
				if(!gl_scale){
					if(fbuf->my_mst_info->my_mst == gl_my_rank)
						clean_iodb(fbuf->my_mst_info->iodb, fbuf->nvars);

					delete_ioreqs(fbuf,0);
				}
				
				DTF_DBG(VERBOSE_DBG_LEVEL, "Finished transfer for %s", fbuf->file_path);
				DTF_DBG(VERBOSE_ERROR_LEVEL, "Time for matching for %s %.4f", fbuf->file_path, MPI_Wtime() - t_start);				
				if(strstr(fbuf->file_path, "hist.d")!=NULL){
					gl_stats.end_mtch_hist = MPI_Wtime()-gl_stats.walltime;
					gl_stats.t_mtch_hist = MPI_Wtime() - t_start;
				}else if(strstr(fbuf->file_path, "anal.d")!=NULL){
					gl_stats.end_mtch_rest = MPI_Wtime()-gl_stats.walltime;
					gl_stats.t_mtch_rest = MPI_Wtime() - t_start;
				}

				fbuf->cur_transfer_epoch++;
			}
			fbuf = fbuf->next;
		}
		progress_comm();
	}
}

int match_ioreqs(file_buffer_t *fbuf)
{

//    double t_all = MPI_Wtime();
//    double t_part1 = MPI_Wtime();
//    double t_part2=MPI_Wtime();
//    double t_part3=MPI_Wtime();
    double t_start;
    int replay = 0, err;
    fname_pattern_t *pat = NULL;
    
	t_start = MPI_Wtime();
	DTF_DBG(VERBOSE_DBG_LEVEL, "time_stamp Match ioreqs for file %s (ncid %d)", fbuf->file_path, fbuf->ncid);
   

	if(strstr(fbuf->file_path, "hist.d")!=NULL)
		gl_stats.st_mtch_hist = MPI_Wtime()-gl_stats.walltime;
	else if(strstr(fbuf->file_path, "anal.d")!=NULL)
		gl_stats.st_mtch_rest = MPI_Wtime()-gl_stats.walltime;


	/*Check if we should do normal matching or replay a previously 
	 * recorded I/O pattern.*/
	pat = find_fname_pattern(fbuf->file_path);
	assert(pat != NULL);
	if(pat->replay_io){
		if(fbuf->writer_id == gl_my_comp_id && pat->wrt_recorded == IO_PATTERN_RECORDED){
			DTF_DBG(VERBOSE_DBG_LEVEL, "1");
			replay = 1;
		}
		else if(fbuf->reader_id == gl_my_comp_id && pat->rdr_recorded == IO_PATTERN_RECORDED){
			DTF_DBG(VERBOSE_DBG_LEVEL, "2");
			replay = 1;
		}
	}

	//TODO this needs to be moved somewhere else
    if(gl_my_comp_id == fbuf->writer_id){ /*Writer*/
		/*First of all make sure all processes have info about the other component*/
		if(!fbuf->cpl_info_shared){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Broadcast info about the other component");
			 if(fbuf->root_writer == gl_my_rank){
				while(fbuf->root_reader == -1)
					progress_comm();
				
			}
			if(!gl_scale){
				err = MPI_Bcast(&(fbuf->cpl_mst_info->nmasters), 1, MPI_INT, 0, fbuf->comm);
				CHECK_MPI(err);
				
				if(fbuf->root_writer != gl_my_rank){
					assert(fbuf->cpl_mst_info->masters== NULL);
					fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
					assert(fbuf->cpl_mst_info->masters != NULL);	
				}
				
				err = MPI_Bcast(fbuf->cpl_mst_info->masters, fbuf->cpl_mst_info->nmasters, MPI_INT, 0, fbuf->comm);
				CHECK_MPI(err);	
			}
			fbuf->root_reader = fbuf->cpl_mst_info->masters[0];
			fbuf->cpl_info_shared = 1;
		 }
	}
	
	if(replay){
		assert(0); //TODO
		DTF_DBG(VERBOSE_DBG_LEVEL, "Will replay I/O instead of matching");
		/*First writer must sync with reader*/
		if(gl_my_comp_id == fbuf->writer_id) 
			 replay_io(pat, fbuf->file_path, fbuf->cur_transfer_epoch);
			
		while(!fbuf->done_matching_flag){
			progress_comm();
		}
			
	} else {
		
		/*If a writer process doesn't have any io requests, it still has to
		  wait for the master process to let it complete.
		  If a reader process does not have any read requests,
		  it notifies the master that it completed matching and returns.*/
		DTF_DBG(VERBOSE_DBG_LEVEL, "Total %d rreqs and %d wreqs", fbuf->rreq_cnt, fbuf->wreq_cnt);
		if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
			send_ioreqs_by_var(fbuf);
		else //if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
			send_ioreqs_by_block(fbuf);

	//    t_part1 = MPI_Wtime() - t_part1;
	//    t_part2 = MPI_Wtime();

		DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching phase");
		
		while(!fbuf->done_matching_flag){					
			progress_comm();
			if( (fbuf->writer_id == gl_my_comp_id) && (fbuf->my_mst_info->my_mst == gl_my_rank)  )
				do_matching(fbuf);
		}
	}
	
	/*Reset flags*/
	if (fbuf->writer_id == gl_my_comp_id)
        /*Reset flag for future matchings*/
        fbuf->my_mst_info->iodb->updated_flag = 0;
    if(fbuf->my_mst_info->my_mst == gl_my_rank)
		fbuf->my_mst_info->nread_completed = 0;
	fbuf->done_matching_flag = 0;
	fbuf->is_transfering = 0;
	
	if(!gl_scale){
		if(fbuf->my_mst_info->my_mst == gl_my_rank)
			clean_iodb(fbuf->my_mst_info->iodb, fbuf->nvars);

		delete_ioreqs(fbuf,0);
	}
	
	/*Synchronize finishing the matching*/
    //~ if(gl_my_comp_id == fbuf->writer_id){
		//~ MPI_Barrier(fbuf->comm);
		//~ if(fbuf->root_writer == gl_my_rank){
			//~ err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, READ_DONE_CONFIRM_TAG, gl_comps[fbuf->reader_id].comm);
			//~ CHECK_MPI(err);
		//~ }
	//~ } else { /*reader*/
		//~ if(fbuf->root_reader == gl_my_rank){
			//~ char filename[MAX_FILE_NAME];
			
			//~ err = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, READ_DONE_CONFIRM_TAG, gl_comps[fbuf->writer_id].comm, MPI_STATUS_IGNORE);
			//~ CHECK_MPI(err);
			//~ assert(strcmp(filename, fbuf->file_path) == 0);
		//~ }
		//~ MPI_Barrier(fbuf->comm);
	//~ }
	
//    t_part2 = MPI_Wtime() - t_part2;
//    t_part3 = MPI_Wtime();

    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);

	//if(gl_scale)
	DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp Stat: Time for matching %.4f", MPI_Wtime() - t_start);
	{
		//~ int rank, avg=0, my;
		//~ MPI_Comm_rank(fbuf->comm, &rank);
		
		if(strstr(fbuf->file_path, "hist.d")!=NULL){
			gl_stats.end_mtch_hist = MPI_Wtime()-gl_stats.walltime;
			gl_stats.t_mtch_hist = MPI_Wtime() - t_start;
			//~ my = gl_stats.t_mtch_hist;
		}else if(strstr(fbuf->file_path, "anal.d")!=NULL){
			gl_stats.end_mtch_rest = MPI_Wtime()-gl_stats.walltime;
			gl_stats.t_mtch_rest = MPI_Wtime() - t_start;
			//~ my = gl_stats.t_mtch_rest;
		}
		//~ err = MPI_Reduce(&my, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, fbuf->comm);
		//~ CHECK_MPI(err);
		//~ if(strstr(fbuf->file_path, "hist.d")!=NULL)
			//~ gl_stats.t_mtch_hist = avg;
		//~ else if(strstr(fbuf->file_path, "anal.d")!=NULL)
			//~ gl_stats.t_mtch_rest = avg;
	}

	fbuf->cur_transfer_epoch++;

    return 0;
}


/*writer->reader || writer->writer*/
void send_data(file_buffer_t *fbuf, void* buf, int bufsz)
{
    int var_id, rdr_rank, err, i;
    io_req_t *ioreq = NULL;
    dtf_var_t *var = NULL;
    size_t rofft = 0, sofft = 0;
    unsigned char *sbuf = NULL;
    size_t sbufsz = 0;
    int def_el_sz;
    unsigned char *rbuf = (unsigned char*)buf;
    MPI_Offset *start, *count;
    int nblocks_written = 0;
    MPI_Offset *new_count = NULL, *new_start = NULL, *tmp;
    MPI_Offset nelems, fit_nelems;
    int min_bl_ndims;
    size_t min_bl_sz;
   
    rdr_rank = (int)(*(MPI_Offset*)(rbuf+rofft));
    rofft += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Sending data to rank %d", rdr_rank);

    sbuf = (unsigned char*)gl_msg_buf; //dtf_malloc(gl_conf.data_msg_size_limit);//
    assert(sbuf != NULL);
    sbufsz = (int)gl_conf.data_msg_size_limit;

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: use data buf of sz %d", (int)sbufsz);

    while(rofft != (size_t)bufsz){

        if(sofft == 0){
            memcpy(sbuf, fbuf->file_path, MAX_FILE_NAME);
            sofft = MAX_FILE_NAME;
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
        ioreq = fbuf->vars[var_id]->ioreqs;
        while(ioreq != NULL){
            if(ioreq->rw_flag == DTF_WRITE){
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
            max_mem = sbufsz - (MAX_FILE_NAME ) -
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
                memcpy(sbuf, fbuf->file_path, MAX_FILE_NAME);
                sofft = MAX_FILE_NAME;
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
                assert(sofft != MAX_FILE_NAME );
                DTF_DBG(VERBOSE_DBG_LEVEL, "Send msg to %d, size %d at %.3f", rdr_rank,(int)sofft, MPI_Wtime() - gl_stats.walltime);
                /*Send this message*/
                MPI_Request req;
                int dsz = (int)sofft;
                double t_start_comm = MPI_Wtime();
                err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm, &req);
                CHECK_MPI(err);
                err = MPI_Wait(&req, MPI_STATUS_IGNORE);
                CHECK_MPI(err);

                gl_stats.t_comm += MPI_Wtime() - t_start_comm;
                gl_stats.ndata_msg_sent++;
                gl_stats.data_msg_sz += sofft;

                sofft = 0;
            } else {
                int cnt_mismatch = 0;
                int type_mismatch = 0;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Will copy subblock (strt->cnt):");
                for(i=0; i < var->ndims; i++)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld \t orig: \t %lld\t --> %lld)", new_start[i], new_count[i], start[i], count[i]);

                /*Copy the block to send buffer*/
                if(var->dtype != ioreq->dtype){
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

                    if(0==0 && cnt_mismatch==0){//it's a continious block of memory
                        int req_el_sz;
                        MPI_Offset start_cpy_offt;
                        MPI_Type_size(ioreq->dtype, &req_el_sz);
                        start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, new_start) * req_el_sz;
                        if(type_mismatch){
                            convertcpy(ioreq->dtype, var->dtype, (unsigned char*)ioreq->user_buf+start_cpy_offt,
                                       (void*)(sbuf+sofft), (int)fit_nelems);
                        }else{
                            memcpy(sbuf+sofft, (unsigned char*)ioreq->user_buf+start_cpy_offt, fit_nelems*def_el_sz);
                        }
                        if(gl_conf.do_checksum){
							double chsum = compute_checksum(sbuf+sofft, var->ndims, new_count, var->dtype);
							DTF_DBG(VERBOSE_DBG_LEVEL, "chsum contin for req %d is %.4f", ioreq->id, chsum);
						}
                    } else {
						get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                           ioreq->count, new_start, new_count, sbuf+sofft, DTF_READ, type_mismatch);
                        if(gl_conf.do_checksum){
							double chsum = compute_checksum(sbuf+sofft, var->ndims, new_count, var->dtype);
							DTF_DBG(VERBOSE_DBG_LEVEL, "chsum getput for req %d is %.4f", ioreq->id, chsum);
						}
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
		err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_comps[fbuf->reader_id].comm, &req);
		CHECK_MPI(err); 
        err = MPI_Wait(&req, MPI_STATUS_IGNORE);
        CHECK_MPI(err);
        gl_stats.t_comm += MPI_Wtime() - t_start_comm;
        gl_stats.ndata_msg_sent++;
        gl_stats.data_msg_sz += sofft;
    }

	/*Record the pattern if needed*/
	record_io_pat(fbuf->file_path, rdr_rank, buf, bufsz, fbuf->cur_transfer_epoch);

    free(new_start);
    free(new_count);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished sending the data");
   // dtf_free(sbuf, sbufsz); //comment out since use gl_msg_buf
}

/*writer->reader*/
static void recv_data_rdr(file_buffer_t *fbuf, void* buf, int bufsz)
{
    int var_id, i, nelems;
    int def_el_sz, req_el_sz;
    dtf_var_t *var = NULL;
    io_req_t *ioreq = NULL;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
    int nblocks_read = 0;
    int type_mismatch;

    MPI_Offset *start, *count;
    void *data;

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

        /*Find the ioreq that has info about this block*/
        ioreq = fbuf->vars[var_id]->ioreqs;
        while(ioreq != NULL){
            if(ioreq->rw_flag == DTF_READ){
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
        assert(nelems <= ioreq->user_buf_sz - ioreq->get_sz);

        DTF_DBG(VERBOSE_DBG_LEVEL, "Will get %d elems for var %d", nelems, var->id);
        if(var->ndims == 0){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Copy scalar variable");
            get_put_data(var, ioreq->dtype, ioreq->user_buf, NULL,
                            NULL, NULL, NULL, data, DTF_WRITE, type_mismatch);
        } else if(cnt_mismatch==0 ){ //continious block of memory
            MPI_Offset start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, start) * req_el_sz;
            if(type_mismatch){
                convertcpy(var->dtype, ioreq->dtype, data, (unsigned char*)ioreq->user_buf + start_cpy_offt, nelems);
            } else
                memcpy((unsigned char*)ioreq->user_buf + start_cpy_offt, data, nelems*def_el_sz);
			//~ if(gl_conf.do_checksum){
				//~ double chsum = compute_checksum((unsigned char*)ioreq->user_buf + start_cpy_offt, var->ndims, count, ioreq->dtype);
				//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "chsum contin for req %d is %.4f", ioreq->id, chsum);
			//~ }
        } else{
			/*Copy from rbuf to user buffer*/
			get_put_data(var, ioreq->dtype, ioreq->user_buf, ioreq->start,
                                       ioreq->count, start, count, data, DTF_WRITE, type_mismatch);
        }
        gl_stats.nbl++;
        offt += nelems*def_el_sz +(nelems*def_el_sz)%sizeof(MPI_Offset);
        ioreq->get_sz += (MPI_Offset)(nelems*req_el_sz);


        DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, var %d, Got %d (expect %d)", ioreq->id, var_id, (int)ioreq->get_sz, (int)ioreq->user_buf_sz);
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
                DTF_DBG(VERBOSE_DBG_LEVEL, "chsum for req %d is %.4f", ioreq->id, chsum);
            }
            //delete this ioreq
            delete_ioreq(fbuf, var_id, &ioreq);

            if(fbuf->rreq_cnt == 0){
				int err;
				
                DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE:Completed all rreqs for file %s", fbuf->file_path);

				dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, READ_DONE_TAG);
				err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
                          READ_DONE_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
				CHECK_MPI(err);

				//Reader ranks except for root complete the send immediately
				if(gl_my_rank != fbuf->my_mst_info->my_mst){
					err = MPI_Wait(&(msg->req), MPI_STATUS_IGNORE);
					CHECK_MPI(err);
					dtf_free(msg, sizeof(dtf_msg_t));
					fbuf->done_matching_flag = 1;
				} else
					ENQUEUE_ITEM(msg, gl_out_msg_q);
            }
        }
    } //while(bufsz)

    DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: time to extract the data: %.3f (%d blocks)", MPI_Wtime() - t_begin, nblocks_read);
}

/*Send the file header and info about vars to the reader when the writer finishes the def mode*/
void send_file_info(file_buffer_t *fbuf, int reader_root)
{
    void *sbuf;
    MPI_Offset sbuf_sz, err;
	int flag;
	
    DTF_DBG(VERBOSE_DBG_LEVEL, "Will send file info to reader %d", reader_root);
    pack_file_info(fbuf, &sbuf_sz, &sbuf);
    assert(sbuf_sz > 0);
    double t_start = MPI_Wtime();
    dtf_msg_t *msg = new_dtf_msg(sbuf, sbuf_sz, DTF_UNDEFINED, FILE_INFO_TAG);
    err = MPI_Isend(sbuf, (int)sbuf_sz, MPI_BYTE, reader_root, FILE_INFO_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
    CHECK_MPI(err);
    err = MPI_Test(&(msg->req), &flag, MPI_STATUS_IGNORE);
	CHECK_MPI(err);
	ENQUEUE_ITEM(msg, gl_out_msg_q);
    gl_stats.t_comm += MPI_Wtime() - t_start;
}


static void notify_masters(file_buffer_t *fbuf, void *msg, int msgsz, int msgtag)
{
    int i, err;
    MPI_Request *reqs = dtf_malloc(fbuf->my_mst_info->nmasters * sizeof(MPI_Request));
    assert(reqs != NULL);
    
    assert(gl_my_rank == fbuf->root_writer);
    
    for(i = 0; i < fbuf->my_mst_info->nmasters; i++) {
        if(fbuf->my_mst_info->masters[i] == gl_my_rank){
			reqs[i] = MPI_REQUEST_NULL;
            continue;
		}
		if(msgsz == 0)
			err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->masters[i], msgtag, gl_comps[gl_my_comp_id].comm, &reqs[i]);
		else 
			err = MPI_Isend(msg, msgsz, MPI_BYTE, fbuf->my_mst_info->masters[i], msgtag, gl_comps[gl_my_comp_id].comm, &reqs[i]);
        CHECK_MPI(err);
    }
    err = MPI_Waitall(fbuf->my_mst_info->nmasters, reqs, MPI_STATUSES_IGNORE);
    CHECK_MPI(err);
    
    dtf_free(reqs, fbuf->my_mst_info->nmasters);
}

/*Function executed by master writers to notify other writers
  about something defined by mpitag*/
static void notify_workgroup(file_buffer_t *fbuf,  void *msg, int msgsz, int msgtag)
{
    int i, err;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Mst %d will notify workgroup (msgtag %d) for %s", gl_my_rank, msgtag, fbuf->file_path);
    assert(fbuf->my_mst_info->my_mst == gl_my_rank);
    MPI_Request *reqs = dtf_malloc((fbuf->my_mst_info->my_wg_sz - 1)*sizeof(MPI_Request));
    assert(reqs != NULL);
    for(i = 0; i < fbuf->my_mst_info->my_wg_sz - 1; i++){
        if(msgsz == 0)
			err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_wg[i], msgtag, gl_comps[gl_my_comp_id].comm, &reqs[i]);
		else 
			err = MPI_Isend(msg, msgsz, MPI_BYTE, fbuf->my_mst_info->my_wg[i], msgtag, gl_comps[gl_my_comp_id].comm, &reqs[i]);
        CHECK_MPI(err);
        CHECK_MPI(err);
    }
    err = MPI_Waitall(fbuf->my_mst_info->my_wg_sz - 1, reqs, MPI_STATUSES_IGNORE);
    CHECK_MPI(err);
    
    dtf_free(reqs, (fbuf->my_mst_info->my_wg_sz - 1)*sizeof(MPI_Request));
}

static void print_recv_msg(int tag, int src)
{
    switch(tag){
        case FILE_READY_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_READY_TAG from %d", src);
            break;
        case IO_DATA_REQ_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_DATA_REQ_TAG from %d", src);
            break;
        case READ_DONE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag READ_DONE_TAG from %d", src);
            break;
        case READ_DONE_CONFIRM_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag READ_DONE_CONFIRM_TAG from %d", src);
            break;
        case FILE_INFO_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_TAG from %d", src);
            break;
        case FILE_INFO_REQ_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_REQ_TAG from %d", src);
            break;
        case MATCH_DONE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag MATCH_DONE_TAG from %d", src);
            break;
        case IO_REQS_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_REQS_TAG from %d", src);
            break;
        case IO_DATA_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_DATA_TAG from %d", src);
            break;
        case COMP_SYNC_TAG:
			DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag COMP_SYNC_TAG from %d", src);
			break;
        default:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag unknown %d from %d", tag, src);
            assert(0);
    }
}

int parce_msg(int comp, int src, int tag, void *rbuf, int bufsz, int is_queued)
{
	MPI_Status status;
	file_buffer_t *fbuf;
	size_t offt = 0;
	char filename[MAX_FILE_NAME];
	int ret = 1;
	
	memcpy(filename, rbuf, MAX_FILE_NAME);
	offt += MAX_FILE_NAME ;
	fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
	if(fbuf == NULL) goto fn_exit;
	if(!fbuf->is_defined) goto fn_exit;
	
	switch(tag){
			case FILE_INFO_REQ_TAG:
				assert(fbuf->root_writer == gl_my_rank);
				if(gl_scale){
					fbuf->root_reader = src;
					assert(fbuf->root_reader == fbuf->cpl_mst_info->masters[0]);
				} else {
					fbuf->root_reader = src;
					//parse mst info of the other component
					assert(fbuf->cpl_mst_info->nmasters == 0);
					memcpy(&(fbuf->cpl_mst_info->nmasters), (unsigned char*)rbuf+offt, sizeof(int));
					assert(fbuf->cpl_mst_info->nmasters > 0);
					offt+=sizeof(int);
					fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
					assert(fbuf->cpl_mst_info->masters != NULL);
					memcpy(fbuf->cpl_mst_info->masters, (unsigned char*)rbuf+offt, fbuf->cpl_mst_info->nmasters*sizeof(int));
					assert(fbuf->root_reader == fbuf->cpl_mst_info->masters[0]);
				}
				
				DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer for file %s, process the file info request, root reader %d",
						fbuf->file_path, fbuf->cpl_mst_info->masters[0] );
			
				if(fbuf->iomode == DTF_IO_MODE_MEMORY) 
					send_file_info(fbuf, fbuf->root_reader);
				break;
			case IO_REQS_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL, "Received reqs from %d, comp %d (my comp %d), bufsz %d", src, comp,gl_my_comp_id, (int)bufsz);
				parse_ioreqs(fbuf, (unsigned char*)rbuf+offt, bufsz - offt, src, gl_comps[comp].comm);
				break;
			case IO_DATA_REQ_TAG:
				send_data(fbuf, (unsigned char*)rbuf+offt, bufsz - offt);
				break;
			case IO_DATA_TAG:
				recv_data_rdr(fbuf, (unsigned char*)rbuf + offt, bufsz - offt);
				break;
			case READ_DONE_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL, "Recv read done for file %s from %d in comp %d (tot %d)", fbuf->file_path,
						src, comp, fbuf->my_mst_info->nread_completed);
				
				fbuf->my_mst_info->nread_completed++;

				if(gl_my_comp_id == fbuf->writer_id){
					assert(comp == fbuf->reader_id);
					assert(gl_my_rank == fbuf->root_writer);
					if(fbuf->my_mst_info->nread_completed == fbuf->cpl_mst_info->nmasters){
						notify_masters(fbuf, NULL, 0, MATCH_DONE_TAG);
						notify_workgroup(fbuf, NULL, 0, MATCH_DONE_TAG);
						fbuf->done_matching_flag = 1;
						DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE, Done matching at %.3f", MPI_Wtime()-gl_stats.walltime);
					}
					
					
				} else { /*reader comp*/
					
					assert(fbuf->my_mst_info->my_mst == gl_my_rank);
					
					if(fbuf->my_mst_info->nread_completed == fbuf->my_mst_info->my_wg_sz){
							int flag=0;
							DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer root that my workgroup completed reading");
							dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, READ_DONE_TAG);
							int err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, READ_DONE_TAG, 
												gl_comps[fbuf->writer_id].comm, &(msg->req));
							CHECK_MPI(err);
							err = MPI_Test(&(msg->req), &flag, MPI_STATUS_IGNORE);
							CHECK_MPI(err);
							ENQUEUE_ITEM(msg, gl_out_msg_q);
							fbuf->done_matching_flag = 1;
							DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
					}
				}
				break;
			case MATCH_DONE_TAG:
				if(fbuf->my_mst_info->my_mst == gl_my_rank)
					notify_workgroup(fbuf, NULL, 0, MATCH_DONE_TAG);
				DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
				fbuf->done_matching_flag = 1;
				break;
			case COMP_SYNC_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL, "Sync tag for file %s", (char*)rbuf);
				fbuf->sync_comp_flag = 1;
				break;
			case FILE_READY_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL,   "Process FILE_READY notif for %s", (char*)rbuf);
				fbuf->is_ready = 1;
				break;
		   default:
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: unknown tag %d", status.MPI_TAG);
				assert(0);
		}
	dtf_free(rbuf, bufsz);	
	return ret;
	
fn_exit:
	if(!is_queued){
		DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: no record for file %s yet, queue the request to process later", filename); 
		dtf_msg_t *msg = new_dtf_msg(rbuf, bufsz, src, tag);
		ENQUEUE_ITEM(msg, gl_comps[comp].in_msg_q);
	}
	return 0;
}

void progress_transfer()
{
	file_buffer_t *fbuf;
	
	progress_comm();
	fbuf = gl_filebuf_list;
	while(fbuf != NULL){
		if( (fbuf->iomode == DTF_IO_MODE_MEMORY) && 
			!fbuf->ignore_io && 
			(fbuf->writer_id == gl_my_comp_id) && 
			(fbuf->my_mst_info->my_mst == gl_my_rank)  )
				do_matching(fbuf);
				
		fbuf = fbuf->next;
	}
}

void progress_comm()
{
    MPI_Status status;
    int comp, flag, err, bufsz, tag, src;
    gl_stats.nprogress_call++;
    double t_st, t_start_comm;
    void *rbuf;
	t_st = MPI_Wtime();
	progress_recv_queue();
    progress_send_queue();

    for(comp = 0; comp < gl_ncomp; comp++){
        if( gl_comps[comp].comm == MPI_COMM_NULL || comp == gl_my_comp_id){
            continue;
        }

        while(1){
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                gl_stats.idle_time += MPI_Wtime() - t_st;
                break;
            } 
            tag = status.MPI_TAG;
            src = status.MPI_SOURCE;
			MPI_Get_count(&status, MPI_BYTE, &bufsz);
			
			rbuf = dtf_malloc(bufsz);
			assert(rbuf != NULL);
            print_recv_msg(tag, src);
            t_start_comm = MPI_Wtime();
			err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, tag, gl_comps[comp].comm, &status);
			gl_stats.t_comm += MPI_Wtime() - t_start_comm;
			CHECK_MPI(err);
			parce_msg(comp, src, tag, rbuf, bufsz, 0);
        }
    }
    //check messages in my component
    comp = gl_my_comp_id;
    while(1){
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                gl_stats.idle_time += MPI_Wtime() - t_st;
                break;
            }         
			tag = status.MPI_TAG;
			src = status.MPI_SOURCE;
			MPI_Get_count(&status, MPI_BYTE, &bufsz);
			
			rbuf = dtf_malloc(bufsz);
			assert(rbuf != NULL);
            print_recv_msg(tag, src);
            t_start_comm = MPI_Wtime();
			err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, tag, gl_comps[comp].comm, &status);
			gl_stats.t_comm += MPI_Wtime() - t_start_comm;
			CHECK_MPI(err);
			parce_msg(comp, src, tag, rbuf, bufsz, 0);
	}
    gl_stats.t_progr_comm += MPI_Wtime() - t_st;
}



void log_ioreq(file_buffer_t *fbuf,
			  int varid, int ndims,
			  const MPI_Offset *start,
			  const MPI_Offset *count,
			  MPI_Datatype dtype,
			  void *buf,
			  int rw_flag)
{
    int el_sz;
	io_req_log_t *req = dtf_malloc(sizeof(io_req_log_t));
	assert(req != NULL);
	req->var_id = varid;
	req->ndims = ndims;
	req->next = NULL;
	req->id = fbuf->rreq_cnt+fbuf->wreq_cnt;
	req->rw_flag = rw_flag;
	req->dtype = dtype;

	if(rw_flag == DTF_READ)
		fbuf->rreq_cnt++;
	else
		fbuf->wreq_cnt++;

    MPI_Type_size(dtype, &el_sz);

	if(ndims > 0){
        int i;
        req->user_buf_sz = count[0];
        for(i=1;i<ndims;i++)
            req->user_buf_sz *= count[i];

        req->start = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        assert(req->start != NULL);
        memcpy((void*)req->start, (void*)start, sizeof(MPI_Offset)*ndims);
        req->count = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        assert(req->count != NULL);
        memcpy((void*)req->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else{
        req->start = NULL;
        req->count = NULL;
        req->user_buf_sz = el_sz;
    }
	/*buffering*/
	if(gl_conf.buffered_req_match && (rw_flag == DTF_WRITE)){
        req->user_buf = dtf_malloc((size_t)req->user_buf_sz);
        assert(req->user_buf != NULL);
        memcpy(req->user_buf, buf, (size_t)req->user_buf_sz);
    } else
        req->user_buf = buf;
    /*checksum*/
    //~ if( (rw_flag == DTF_WRITE) && gl_conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
        //~ req->checksum = compute_checksum(buf, ndims, count, dtype);
        //~ DTF_DBG(VERBOSE_DBG_LEVEL, "checksum %.4f", req->checksum);
    //~ } else
        //~ req->checksum = 0;

	/*Enqueue the request to the head*/
	if(fbuf->ioreq_log == NULL)
		fbuf->ioreq_log = req;
	else{
		req->next = fbuf->ioreq_log;
		fbuf->ioreq_log = req;
	}
}
