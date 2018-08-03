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

static void delete_ioreq(file_buffer_t *fbuf, int varid, io_req_t **ioreq)
{
	DTF_DBG(VERBOSE_ALL_LEVEL, "Delete req %lu, cur wreqs %lu, rreqs %lu", (*ioreq)->id,
			fbuf->wreq_cnt, fbuf->rreq_cnt);

    dtf_var_t *var = fbuf->vars[varid];

    if( (*ioreq)->is_buffered)
        dtf_free((*ioreq)->user_buf, (size_t) (*ioreq)->req_data_sz);

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
    if((*ioreq)->derived_params != NULL){
        dtf_free((*ioreq)->derived_params->orig_array_size, var->ndims*sizeof(MPI_Offset));
        dtf_free((*ioreq)->derived_params->orig_start, var->ndims*sizeof(MPI_Offset));
        dtf_free((*ioreq)->derived_params, sizeof(dtype_params_t));
    } 
    dtf_free((*ioreq), sizeof(io_req_t));
}

//~ static void invalidate_old_ioreqs(file_buffer_t *fbuf)
//~ {
	//~ io_req_t *ioreq, *tmp;
    //~ int varid;
    //~ dtf_var_t *var;
    
	//~ for(varid=0; varid < fbuf->nvars; varid++){
		//~ var = fbuf->vars[varid];
		//~ ioreq = var->ioreqs;
		//~ while(ioreq != NULL){
			
			//~ if(ioreq->sent_flag) break;
			
			//~ if((fbuf->writer_id == gl_proc.my_comp) && (ioreq->rw_flag == DTF_READ)){
				//~ DTF_DBG(VERBOSE_ALL_LEVEL, "Invalidate ioreq %d", ioreq->id);
				//~ assert(!ioreq->sent_flag);
				//~ tmp = ioreq->next;
				//~ delete_ioreq(fbuf, varid, &ioreq);
				//~ ioreq = tmp;
			//~ } else 
				//~ ioreq = ioreq->next;
		//~ }
	//~ }
//~ }

void delete_ioreqs(file_buffer_t *fbuf)
{
    io_req_t *ioreq, *tmp;
    int varid;
    dtf_var_t *var;
    DTF_DBG(VERBOSE_ALL_LEVEL, "Delete io requests for file %s", fbuf->file_path);
	for(varid=0; varid < fbuf->nvars; varid++){
		var = fbuf->vars[varid];
		ioreq = var->ioreqs;
		while(ioreq != NULL){
			tmp = ioreq->next;
			delete_ioreq(fbuf, varid, &ioreq);
			ioreq = tmp;
		}
	}
	assert(fbuf->rreq_cnt == 0);
	assert(fbuf->wreq_cnt == 0);
	
}

/*match subblocks of data*/
static void do_matching(file_buffer_t *fbuf)
{
    int mlc_ranks = 4;
    int mlc_buf   = 1024;
    int i, j;
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
    sbuf = (unsigned char**)dtf_malloc(mlc_ranks*sizeof(unsigned char*));
    bufsz = (int*)dtf_malloc(mlc_ranks*sizeof(int));
    offt = (size_t*)dtf_malloc(mlc_ranks*sizeof(size_t));

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
    for(j = 0; j < fbuf->cpl_mst_info->comm_sz; j++){
		
		if(fbuf->my_mst_info->iodb->ritems[j] == NULL) continue;
		
		ritem = fbuf->my_mst_info->iodb->ritems[j];

        t_st = MPI_Wtime();

        n_matched_blocks = 0;

        for(i = 0; i < allocd_nwriters; i++){
            sbuf[i] = NULL;
            bufsz[i] = 0;
            offt[i] = 0;
            //writers[i] = -1;
        }

        nwriters = 0;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Matching rreq from rank %d", ritem->global_rank);
        //~ if(ritem->comm == gl_proc.comps[gl_proc.my_comp].comm)
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
				gl_proc.stats_info.idle_time += MPI_Wtime() - t_st;
				t_idle += MPI_Wtime() - t_st;
				gl_proc.stats_info.idle_do_match_time += MPI_Wtime() - t_st;
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
            } else
                matched_count = NULL;

   
            while(nelems_to_match){

                nelems_matched = 0;
				if(var->ndims > 0)
					wblock = rb_find_block(witem->dblocks, rblock->start, rblock->count, var->ndims);
				else
					wblock = (block_t*)witem->dblocks;
				
                if(wblock == NULL){

                    //didn't find
                    DTF_DBG(VERBOSE_ALL_LEVEL, "didnt find block for var %d", var_id);
                    gl_proc.stats_info.idle_time += MPI_Wtime() - t_st;
                    t_idle += MPI_Wtime() - t_st;
                    gl_proc.stats_info.idle_do_match_time += MPI_Wtime() - t_st;
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
                            gl_proc.stats_info.malloc_size += mlc_ranks*sizeof(int);

                            tmp = realloc((void*)offt, (nwriters+mlc_ranks)*sizeof(size_t));
                            assert(tmp != NULL);
                            offt = (size_t*)tmp;
                            gl_proc.stats_info.malloc_size += mlc_ranks*sizeof(size_t);

                            tmp = realloc((void*)bufsz, (nwriters+mlc_ranks)*sizeof(int));
                            assert(tmp != NULL);
                            bufsz = (int*)tmp;
                            gl_proc.stats_info.malloc_size += mlc_ranks*sizeof(int);

                            tmp1 = realloc(sbuf, (nwriters+mlc_ranks)*sizeof(unsigned char*));
                            assert(tmp1 != NULL);
                            sbuf = tmp1;
                            gl_proc.stats_info.malloc_size += mlc_ranks*sizeof(unsigned char*);
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
						*(MPI_Offset*)(sbuf[rank_idx] + offt[rank_idx]) = (MPI_Offset)ritem->global_rank; //to whom writer should send the data
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
                if(writers[i] == gl_proc.myrank){
                    my_idx = i;
                    continue;
                } else {
                    dtf_msg_t *msg = new_dtf_msg(sbuf[i], bufsz[i], DTF_UNDEFINED, IO_DATA_REQ_TAG, 1);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Send data req to wrt %d", writers[i]);
                    err = MPI_Isend((void*)sbuf[i], offt[i], MPI_BYTE, writers[i], IO_DATA_REQ_TAG, gl_proc.comps[gl_proc.my_comp].comm, msg->reqs);
                    CHECK_MPI(err);
                    ENQUEUE_ITEM(msg, gl_proc.comps[gl_proc.my_comp].out_msg_q);
                }
            }

            if(my_idx != -1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Parse data req to myself");
                send_data(fbuf, sbuf[my_idx] + MAX_FILE_NAME, offt[my_idx] - MAX_FILE_NAME);
                dtf_free(sbuf[my_idx], bufsz[my_idx]);
            }
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched rreq for rank %d", ritem->global_rank);
        }

        /*If we matched all chunks for this rank, then delete this ritem*/
        if(ritem->nblocks == 0){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Matched all. Delete ritem of rank %d (left ritems %d). ", ritem->global_rank,  (int)(fbuf->my_mst_info->iodb->nritems - 1));
            dtf_free(ritem, sizeof(read_db_item_t));
            fbuf->my_mst_info->iodb->ritems[j] = NULL;
            fbuf->my_mst_info->iodb->nritems--;
        } else {
           DTF_DBG(VERBOSE_DBG_LEVEL, "Haven't matched all blocks for rank %d, %d left", ritem->global_rank, (int)ritem->nblocks);
        }
    }
    if(nwriters == 0)
        t_st = MPI_Wtime(); //reset
    /*dealloc stuff*/
    dtf_free(writers, allocd_nwriters*sizeof(int));
    dtf_free(offt, allocd_nwriters*sizeof(size_t));
    dtf_free(bufsz, allocd_nwriters*sizeof(int));
    dtf_free(sbuf, allocd_nwriters*sizeof(unsigned char*));
    gl_proc.stats_info.ndb_match++;

    if(nwriters == 0)
        gl_proc.stats_info.idle_time += MPI_Wtime() - t_st;

    assert( MPI_Wtime() - t_start >= t_idle);
    gl_proc.stats_info.master_time = MPI_Wtime() - t_start - t_idle;  //useful work
	gl_proc.stats_info.t_do_match += MPI_Wtime() - t_start;

    DTF_DBG(VERBOSE_DBG_LEVEL, "after matching: %d ritems", (int)fbuf->my_mst_info->iodb->nritems);
}


static void parse_ioreqs(file_buffer_t *fbuf, void *buf, int bufsz, int global_rank, MPI_Comm comm)
{
    int var_id, rw_flag, i;
    dtf_var_t *var = NULL;
    size_t offt = 0;
    unsigned char *chbuf = (unsigned char*)buf;
	read_db_item_t *ritem = NULL;
	int file_rank;
    double t_st = MPI_Wtime();
    
   // offt += MAX_FILE_NAME;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Start parsing reqs for file %s", fbuf->file_path);
    if(comm == gl_proc.comps[fbuf->reader_id].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Reqs are from reader");
    else {
        assert(comm == gl_proc.comps[fbuf->writer_id].comm);
        DTF_DBG(VERBOSE_DBG_LEVEL, "Req are from writer");
    }
	
	if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
	   ((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Discard ioreqs as the other component has started finalizing.");
			return;
	}
	if(fbuf->done_matching_flag){
		DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: received io reqs for file %s after matching already finished, discard ioreqs", fbuf->file_path);
		return;
	}
	
	file_rank = (int)(*(MPI_Offset*)(chbuf+offt));
	offt += sizeof(MPI_Offset);
	
    while(offt != (size_t)bufsz){
        rw_flag = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        var_id = (int)(*(MPI_Offset*)(chbuf+offt));
        offt += sizeof(MPI_Offset);
        gl_proc.stats_info.iodb_nioreqs++;
        assert(var_id < fbuf->nvars);
		var = fbuf->vars[var_id];

        if(rw_flag == DTF_READ){
			//assert(comm == gl_proc.comps[fbuf->reader_id].comm);
			if(fbuf->my_mst_info->iodb->ritems == NULL){
				assert(fbuf->cpl_mst_info->comm_sz > 0);

				fbuf->my_mst_info->iodb->ritems = dtf_malloc(fbuf->cpl_mst_info->comm_sz * sizeof(read_db_item_t*));
				for(i = 0; i < fbuf->cpl_mst_info->comm_sz; i++)
					fbuf->my_mst_info->iodb->ritems[i] = NULL;
			}
			if(file_rank >= fbuf->cpl_mst_info->comm_sz)
				DTF_DBG(VERBOSE_DBG_LEVEL, "comm_sz %d, file_rank %d", fbuf->cpl_mst_info->comm_sz, file_rank);
			assert(file_rank < fbuf->cpl_mst_info->comm_sz);
			ritem = fbuf->my_mst_info->iodb->ritems[file_rank];
			
			if(ritem == NULL){
				fbuf->my_mst_info->iodb->ritems[file_rank] = (read_db_item_t*)dtf_malloc(sizeof(read_db_item_t));
				ritem = fbuf->my_mst_info->iodb->ritems[file_rank];
				ritem->dblocks = NULL;
				ritem->nblocks = 0;
				ritem->global_rank = global_rank;
				fbuf->my_mst_info->iodb->nritems++;
			}
			
            read_dblock_t *dblock = dtf_malloc(sizeof(read_dblock_t));
            dblock->var_id = var_id;
            dblock->ndims = var->ndims;
            dblock->next = NULL;
            dblock->prev = NULL;
            dblock->start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
            dblock->count = dtf_malloc(var->ndims*sizeof(MPI_Offset));
            memcpy(dblock->start, chbuf+offt, var->ndims*sizeof(MPI_Offset));
            offt += var->ndims*sizeof(MPI_Offset);
            memcpy(dblock->count, chbuf+offt, var->ndims*sizeof(MPI_Offset));
            offt += var->ndims*sizeof(MPI_Offset);
            
            DTF_DBG(VERBOSE_ALL_LEVEL, "Added dblock:");
			for(i = 0; i < var->ndims; i++)
				DTF_DBG(VERBOSE_ALL_LEVEL, "%lld  ->  %lld", dblock->start[i], dblock->count[i]);
			
            /*add to list*/
            if(ritem->dblocks == NULL){
                ritem->dblocks = dblock;
            } else {
                dblock->next = ritem->dblocks;
                ritem->dblocks->prev = dblock;
                ritem->dblocks = dblock;
            }
            ritem->nblocks++;
            DTF_DBG(VERBOSE_DBG_LEVEL, "ritems %d, cur item (r %d) %lld blocks",(int)fbuf->my_mst_info->iodb->nritems, global_rank, ritem->nblocks);
        } else { /*DTF_WRITE*/
			int i;
            write_db_item_t *witem;
            /*Allow write requests only from the writer*/
            //assert(comm == gl_proc.comps[fbuf->writer_id].comm);
            
            if(fbuf->my_mst_info->iodb->witems == NULL){
				fbuf->my_mst_info->iodb->witems = dtf_malloc(fbuf->nvars*sizeof(write_db_item_t*));
				for(i = 0; i < fbuf->nvars; i++)
					fbuf->my_mst_info->iodb->witems[i] = NULL;
			}
            if(fbuf->my_mst_info->iodb->witems[var_id] == NULL){
                fbuf->my_mst_info->iodb->witems[var_id] = (write_db_item_t*)dtf_malloc(sizeof(write_db_item_t));
                witem = fbuf->my_mst_info->iodb->witems[var_id];

                if(var->ndims > 0)     
                    witem->dblocks = (void*)RBTreeCreateBlocks(rb_key_cmp, NullFunction, rb_destroy_node_info, rb_print_key, rb_print_info, 0);
                else{
                    witem->dblocks = dtf_malloc(sizeof(block_t));
				}
                witem->nblocks = 0;
                witem->ndims = var->ndims;
            } else
				witem = fbuf->my_mst_info->iodb->witems[var_id];
				
            if(var->ndims > 0){
				insert_info *info = dtf_malloc(sizeof(insert_info));
				info->ndims = var->ndims;

				info->blck = dtf_malloc(sizeof(block_t));
				info->cur_dim = 0;
				info->blck->rank = global_rank;
				info->blck->start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                info->blck->count = dtf_malloc(var->ndims*sizeof(MPI_Offset)); 
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
                ((block_t*)witem->dblocks)->rank = global_rank;
                DTF_DBG(VERBOSE_DBG_LEVEL, "Inserted scalar var");
            }
            witem->nblocks++;
        }
    }
    assert(offt == (size_t)bufsz);
	fbuf->my_mst_info->iodb->updated_flag = 1;
    gl_proc.stats_info.parse_ioreq_time += MPI_Wtime() - t_st;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished parsing reqs. (mem %lu)", gl_proc.stats_info.malloc_size);
}

io_req_t *new_ioreq(int id,
                    int ndims,
                    MPI_Datatype dtype,
                    const MPI_Offset *start,
                    const MPI_Offset *count,
                    dtype_params_t *derived_params,
                    void *buf,
                    int rw_flag,
                    int buffered)
{
    io_req_t *ioreq = (io_req_t*)dtf_malloc(sizeof(io_req_t));
    int el_sz;
	
	MPI_Type_size(dtype, &el_sz);
	
    if(ndims > 0){
        int i;
        ioreq->req_data_sz = count[0];
        for(i=1;i<ndims;i++)
            ioreq->req_data_sz *= count[i];
        ioreq->req_data_sz *= el_sz;

        ioreq->start = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        memcpy((void*)ioreq->start, (void*)start, sizeof(MPI_Offset)*ndims);
        ioreq->count = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        memcpy((void*)ioreq->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else {
        ioreq->start = NULL;
        ioreq->count = NULL;
        ioreq->req_data_sz = el_sz;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "req %d, user bufsz %d", id, (int)ioreq->req_data_sz);
    
    ioreq->derived_params = derived_params;		
    ioreq->is_buffered = buffered;
	ioreq->next = NULL;
    ioreq->sent_flag = 0;
    ioreq->id = id;
    ioreq->get_sz = 0;
    ioreq->prev = NULL;
    ioreq->rw_flag = rw_flag;
    ioreq->dtype = dtype;
    ioreq->user_buf = buf;
    
    if(buffered && (rw_flag == DTF_WRITE)){
        double t_start = MPI_Wtime();
        ioreq->user_buf = dtf_malloc((size_t)ioreq->req_data_sz);
        if(derived_params == NULL)
			memcpy(ioreq->user_buf, buf, (size_t)ioreq->req_data_sz);
		else{
			//Buffer data contiguously and delete the derived parameters since
			//won't need them anymore	 
			get_put_data(ndims, dtype, dtype, derived_params->orig_array_size, buf, 
					     derived_params->orig_start, ioreq->count, ioreq->user_buf, DTF_READ, 0);
			
			dtf_free(derived_params->orig_array_size, ndims * sizeof(MPI_Offset));
			dtf_free(derived_params->orig_start, ndims * sizeof(MPI_Offset));
			dtf_free(derived_params, sizeof(dtype_params_t));
			ioreq->derived_params = NULL;
        }
        gl_proc.stats_info.accum_dbuff_time += MPI_Wtime() - t_start;
		gl_proc.stats_info.accum_dbuff_sz += (size_t)ioreq->req_data_sz;
    } 

    if( (rw_flag == DTF_WRITE) && gl_proc.conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
        ioreq->checksum = compute_checksum(buf, ndims, count, dtype);
        DTF_DBG(VERBOSE_DBG_LEVEL, "chsum for req %lu is %.4f", ioreq->id, ioreq->checksum);
    } else
        ioreq->checksum = 0;
        
    gl_proc.stats_info.nioreqs++;
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
	int file_rank;
	
    if(fbuf->rreq_cnt == 0 && fbuf->wreq_cnt == 0)
        return;

    if(gl_proc.my_comp == fbuf->writer_id)
        mst_info = fbuf->my_mst_info;
    else
        mst_info = fbuf->cpl_mst_info;
    nmasters = mst_info->nmasters;

	MPI_Comm_rank(fbuf->comm, &file_rank);
	
    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));

    /*Distribute ioreqs between the masters based on var id*/

    //alloc mem
    for(mst = 0; mst < nmasters; mst++){
        bufsz[mst] = 0;
        sbuf[mst] = NULL;
    }

	//~ invalidate_old_ioreqs(fbuf);
	
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

    if(fbuf->reader_id == gl_proc.my_comp && nrreqs == 0){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
		int flag; 
		char *fname = dtf_malloc(MAX_FILE_NAME);
		strcpy(fname, fbuf->file_path);
		dtf_msg_t *msg = new_dtf_msg(fname, MAX_FILE_NAME, DTF_UNDEFINED, READ_DONE_TAG, 1);
		err = MPI_Isend(fname, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
				  READ_DONE_TAG, gl_proc.comps[fbuf->reader_id].comm, msg->reqs);
		CHECK_MPI(err);
		err = MPI_Test(msg->reqs, &flag, MPI_STATUS_IGNORE);
		CHECK_MPI(err);
		ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->reader_id].out_msg_q);
		
		if(gl_proc.myrank != fbuf->root_reader){
			fbuf->done_matching_flag = 1;
		}	
    }
    if(!data_to_send)
        goto fn_exit;   //nothing to send

    for(mst = 0; mst < nmasters; mst++){
        if(bufsz[mst] == 0)
            continue;

        bufsz[mst] += MAX_FILE_NAME +sizeof(MPI_Offset);
        DTF_DBG(VERBOSE_DBG_LEVEL, "bufsz %lu for mst %d", bufsz[mst], mst);
        sbuf[mst] = dtf_malloc(bufsz[mst]);

        memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
        offt[mst] = MAX_FILE_NAME ;
		*(MPI_Offset*)(sbuf[mst] + offt[mst]) = file_rank;
		offt[mst] += sizeof(MPI_Offset);
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
			DTF_DBG(VERBOSE_DBG_LEVEL, "Will send ioreq %lu (varid %d) to mst %d", ioreq->id, varid, mst); 
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
        if( (fbuf->writer_id == gl_proc.my_comp) && (mst_info->masters[mst] == gl_proc.myrank)){
            idx = mst;
        } else {
            int flag;
            MPI_Status status;
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], DTF_UNDEFINED, IO_REQS_TAG, 1);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu), comp id %d (mine %d)", mst_info->masters[mst], offt[mst], bufsz[mst], fbuf->writer_id, gl_proc.my_comp);
            assert(offt[mst] == bufsz[mst]);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, mst_info->masters[mst], IO_REQS_TAG,
                            gl_proc.comps[fbuf->writer_id].comm, msg->reqs);
            CHECK_MPI(err);
            err = MPI_Test(msg->reqs, &flag, &status);
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->writer_id].out_msg_q);
        }
    }
    
    if(idx != -1){
		assert(offt[idx] == bufsz[idx]);
        parse_ioreqs(fbuf,sbuf[idx], offt[idx], gl_proc.myrank, gl_proc.comps[gl_proc.my_comp].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    } 
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
    int file_rank;

    if(fbuf->rreq_cnt == 0 && fbuf->wreq_cnt==0)
        return;

    if(gl_proc.my_comp == fbuf->writer_id)
        mst_info = fbuf->my_mst_info;
    else
        mst_info = fbuf->cpl_mst_info;
        
    MPI_Comm_rank(fbuf->comm, &file_rank);

    nmasters = mst_info->nmasters;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Sending I/O reqs by block");

    sbuf = (unsigned char**)dtf_malloc(nmasters*sizeof(unsigned char*));
    offt = (size_t*)dtf_malloc(nmasters*sizeof(size_t));
    bufsz = (size_t*)dtf_malloc(nmasters*sizeof(size_t));

    /*Distribute ioreqs between the masters based on the first dimension of the var*/

    for(mst = 0; mst < nmasters; mst++){
        bufsz[mst] = 0;
        offt[mst] = 0;
        sbuf[mst] = NULL;
    }

	//~ invalidate_old_ioreqs(fbuf);

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
   if(fbuf->reader_id == gl_proc.my_comp && nrreqs == 0){
		 DTF_DBG(VERBOSE_DBG_LEVEL, "Have no read requests. Notify master that read done");
         int flag = 0;
		char *fname = dtf_malloc(MAX_FILE_NAME);
		strcpy(fname, fbuf->file_path);
		dtf_msg_t *msg = new_dtf_msg(fname, MAX_FILE_NAME, DTF_UNDEFINED, READ_DONE_TAG, 1);
		err = MPI_Isend(fname, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
				  READ_DONE_TAG, gl_proc.comps[fbuf->reader_id].comm, msg->reqs);
		CHECK_MPI(err);
		err = MPI_Test(msg->reqs, &flag, MPI_STATUS_IGNORE);
		CHECK_MPI(err);
		ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->reader_id].out_msg_q);
		
		if(gl_proc.myrank != fbuf->root_reader){
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
			DTF_DBG(VERBOSE_DBG_LEVEL, "Will send ioreq for var %d:", varid);
			for(i = 0; i < var->ndims; i++)
				DTF_DBG(VERBOSE_ALL_LEVEL, "%lld  ->  %lld", ioreq->start[i], ioreq->count[i]);
			/*I/O reqs for scalar vars and 1D vars are distributed in round-robin fashion among
			 * matchers. I/O reqs for higher dimension vars are split by among matchers along 
			 * the slowest changing dimension*/
			if(var->ndims == 0 || var->ndims == 1){
				mst = varid % nmasters;
				if(bufsz[mst] == 0){
						sbuf[mst] = dtf_malloc(MAX_FILE_NAME + mlc_chunk);
						bufsz[mst] = MAX_FILE_NAME + mlc_chunk;

						memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
						offt[mst] = MAX_FILE_NAME ;
						*(MPI_Offset*)(sbuf[mst] + offt[mst]) = file_rank;
						offt[mst] += sizeof(MPI_Offset);
				}
				
				if(var->ndims*sizeof(MPI_Offset)*3+offt[mst] > bufsz[mst]){
					size_t ext_sz = mlc_chunk;
					while(bufsz[mst]+ext_sz < var->ndims*sizeof(MPI_Offset)*3+offt[mst])
						ext_sz += mlc_chunk;
					DTF_DBG(VERBOSE_ALL_LEVEL, "bufsz %lu, ext sz %lu", bufsz[mst], ext_sz );
					
					void *tmp = realloc(sbuf[mst], bufsz[mst] + ext_sz);
					assert(tmp != NULL);
					gl_proc.stats_info.malloc_size += ext_sz;
					bufsz[mst] += ext_sz;
				}


				if(var->ndims == 0){
					/*Only save the rw flag and the var id*/
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
					offt[mst] += sizeof(MPI_Offset);
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)varid;
					offt[mst] += sizeof(MPI_Offset);
				} else {
					/*Store var_id, rw_flag, start[] and count[]*/
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)ioreq->rw_flag;
					offt[mst] += sizeof(MPI_Offset);
					*(MPI_Offset*)(sbuf[mst] + offt[mst]) = (MPI_Offset)varid;
					offt[mst] += sizeof(MPI_Offset);
					memcpy(sbuf[mst]+offt[mst], ioreq->start, var->ndims*sizeof(MPI_Offset));
					offt[mst] += var->ndims*sizeof(MPI_Offset);
					memcpy(sbuf[mst]+offt[mst], ioreq->count, var->ndims*sizeof(MPI_Offset));
					offt[mst] += var->ndims*sizeof(MPI_Offset);
				}	
				block_cnt++;
			} else {

				if(var->shape[0] == DTF_UNLIMITED)
					block_range = DEFAULT_BLOCK_SZ_RANGE;
				else if(gl_proc.conf.iodb_range > 0)
					block_range = gl_proc.conf.iodb_range;
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
					if( (block_range == DEFAULT_BLOCK_SZ_RANGE) || (block_range == gl_proc.conf.iodb_range) )
						mst = (int)( ((ioreq->start[0] + shift)/block_range) % nmasters);
					else
						mst = (int)((ioreq->start[0] + shift)/block_range);

					DTF_DBG(VERBOSE_ALL_LEVEL, "mst %d (%lld)", mst, ioreq->start[0] + shift);
					assert(mst < nmasters);
					if(bufsz[mst] == 0){
						sbuf[mst] = dtf_malloc(mlc_chunk);
						bufsz[mst] = mlc_chunk;

						memcpy(sbuf[mst], fbuf->file_path, MAX_FILE_NAME);
						offt[mst] = MAX_FILE_NAME;
						*(MPI_Offset*)(sbuf[mst] + offt[mst]) = file_rank;
						offt[mst] += sizeof(MPI_Offset);
					}

					if(var->ndims*sizeof(MPI_Offset)*3+offt[mst] > bufsz[mst]){
						size_t ext_sz = mlc_chunk;
						while(bufsz[mst]+ext_sz < var->ndims*sizeof(MPI_Offset)*3+offt[mst])
							ext_sz += mlc_chunk;
						DTF_DBG(VERBOSE_ALL_LEVEL, "bufsz %lu, ext sz %lu", bufsz[mst], ext_sz );
						
						void *tmp = realloc(sbuf[mst], bufsz[mst] + ext_sz);
						assert(tmp != NULL);
						gl_proc.stats_info.malloc_size += ext_sz;
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
        if( (fbuf->writer_id == gl_proc.my_comp) && (mst_info->masters[mst] == gl_proc.myrank)){
            idx = mst;
        } else {
            int flag;
            MPI_Status status;
            dtf_msg_t *msg = new_dtf_msg(sbuf[mst], bufsz[mst], DTF_UNDEFINED, IO_REQS_TAG, 1);
            DTF_DBG(VERBOSE_DBG_LEVEL, "Post send ioreqs req to mst %d (bufsz %lu (allcd %lu), comp id %d (mine %d)", mst_info->masters[mst], 		offt[mst], bufsz[mst], fbuf->writer_id, gl_proc.my_comp);
            err = MPI_Isend((void*)sbuf[mst], (int)offt[mst], MPI_BYTE, mst_info->masters[mst], IO_REQS_TAG,
                            gl_proc.comps[fbuf->writer_id].comm, msg->reqs);
            CHECK_MPI(err);
            err = MPI_Test(msg->reqs, &flag, &status);
            CHECK_MPI(err);
            ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->writer_id].out_msg_q);
        }
    }

    //gl_proc.stats_info.t_comm += MPI_Wtime() - t_start_comm;

    if(idx != -1){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Parse reqs to myself: %d", (int)offt[idx] - MAX_FILE_NAME);
        parse_ioreqs(fbuf, (unsigned char*)sbuf[idx] + MAX_FILE_NAME, offt[idx] - MAX_FILE_NAME, gl_proc.myrank, gl_proc.comps[gl_proc.my_comp].comm);
        dtf_free(sbuf[idx], bufsz[idx]);
    }
    
fn_exit:
    dtf_free(sbuf, nmasters*sizeof(unsigned char*));
    dtf_free(bufsz, nmasters*sizeof(unsigned char*));
    dtf_free(offt, nmasters*sizeof(unsigned char*));
}

void match_ioreqs_all_files()
{
	file_buffer_t *fbuf;
	int file_cnt = 0;
	double t_start = MPI_Wtime();
	
	/*First check for how many files need to complete transfer 
	 * and send any unsent I/O requests*/
	fbuf = gl_proc.filebuf_list;
	while(fbuf != NULL){
		if(fbuf->is_transferring){ 
			DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is in active transfer", fbuf->file_path);
			
			if(!fbuf->done_matching_flag && (((fbuf->writer_id == gl_proc.my_comp) && !gl_proc.comps[fbuf->reader_id].finalized)  ||
				((fbuf->reader_id == gl_proc.my_comp) && !gl_proc.comps[fbuf->writer_id].finalized)) ){
				
				if(gl_proc.conf.iodb_build_mode == IODB_BUILD_VARID)
					send_ioreqs_by_var(fbuf);
				else //if(gl_proc.conf.iodb_build_mode == IODB_BUILD_BLOCK)
					send_ioreqs_by_block(fbuf);
			}
		}
		fbuf = fbuf->next;
	}
	
	while(1){
		
		progress_comm();
		
		file_cnt = 0;
		fbuf = gl_proc.filebuf_list;
		while(fbuf != NULL){
			
			if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
				((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: skipping a data transfer session for file %s as the other component has finalized", fbuf->file_path);
				fbuf->done_matching_flag = 1;
			}
			
			if(fbuf->is_transferring && fbuf->done_matching_flag){
				
				/*Reset flags*/
				if(fbuf->my_mst_info->my_mst == gl_proc.myrank)
					fbuf->my_mst_info->nread_completed = 0;
				fbuf->done_matching_flag = 0;
				fbuf->is_transferring = 0;
				
				DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);
				DTF_DBG(VERBOSE_ERROR_LEVEL, "dtf_time transfer for %s: %.4f", fbuf->file_path, MPI_Wtime() - t_start);

				fbuf->cur_transfer_epoch++;
					
				fname_pattern_t *pat = find_fname_pattern(fbuf->file_path);			
				assert(pat != NULL);
				if(fbuf->session_cnt == pat->num_sessions){
					file_buffer_t *tmp = fbuf->next;
					delete_file_buffer(fbuf);
					fbuf = tmp;
					continue;
				}
				
			}
			
			if(fbuf->is_transferring) file_cnt++;
			
			do_matching(fbuf);
						
			fbuf = fbuf->next;
		}
		if(file_cnt == 0) break;
	}
}

int match_ioreqs(file_buffer_t *fbuf)
{

//    double t_all = MPI_Wtime();
//    double t_part1 = MPI_Wtime();
//    double t_part2=MPI_Wtime();
//    double t_part3=MPI_Wtime();
    double t_start;
    int replay = 0;
    fname_pattern_t *pat = NULL;
    
	t_start = MPI_Wtime();
	DTF_DBG(VERBOSE_DBG_LEVEL, "Match ioreqs for file %s (ncid %d)", fbuf->file_path, fbuf->ncid);

	/*Check if we should do normal matching or replay a previously 
	 * recorded I/O pattern.*/
	pat = find_fname_pattern(fbuf->file_path);
	assert(pat != NULL);
	if(pat->replay_io){
		if(fbuf->writer_id == gl_proc.my_comp && pat->wrt_recorded == IO_PATTERN_RECORDED){
			replay = 1;
		}
		else if(fbuf->reader_id == gl_proc.my_comp && pat->rdr_recorded == IO_PATTERN_RECORDED){
			replay = 1;
		}
	}
	
	if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
		((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: skipping a data transfer session for file %s as the other component has finalized", fbuf->file_path);
		fbuf->done_matching_flag = 1;
		goto done_match;
	}

	if(replay){
		
		while(fbuf->cpl_mst_info->nmasters == 0){
			progress_comm();
			if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
				((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: skipping a data transfer session for file %s as the other component has finalized", fbuf->file_path);
				fbuf->done_matching_flag = 1;
				goto done_match;
			}
		}
		
		DTF_DBG(VERBOSE_DBG_LEVEL, "Will replay I/O instead of matching");
		
		if(gl_proc.my_comp == fbuf->writer_id){ 
			 replay_io_pat(pat, fbuf->file_path, fbuf->cur_transfer_epoch);
			 fbuf->done_matching_flag = 1;
		} else 		
			while(!fbuf->done_matching_flag) progress_comm();
			
	} else {

			/*If a writer process doesn't have any io requests, it still has to
			  wait for the master process to let it complete.
			  If a reader process does not have any read requests,
			  it notifies the master that it completed matching and returns.*/
			DTF_DBG(VERBOSE_DBG_LEVEL, "Total %lu rreqs and %lu wreqs", fbuf->rreq_cnt, fbuf->wreq_cnt);
			if(gl_proc.conf.iodb_build_mode == IODB_BUILD_VARID)
				send_ioreqs_by_var(fbuf);
			else //if(gl_proc.conf.iodb_build_mode == IODB_BUILD_BLOCK)
				send_ioreqs_by_block(fbuf);

			DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching phase");
			
			while(!fbuf->done_matching_flag){					
				progress_comm();
				if( (fbuf->writer_id == gl_proc.my_comp) && (fbuf->my_mst_info->my_mst == gl_proc.myrank)  )
					do_matching(fbuf);
				
				if( ((fbuf->writer_id == gl_proc.my_comp) && gl_proc.comps[fbuf->reader_id].finalized)  ||
					((fbuf->reader_id == gl_proc.my_comp) && gl_proc.comps[fbuf->writer_id].finalized) ){
					DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: skipping a data transfer session for file %s as the other component has finalized", fbuf->file_path);
					fbuf->done_matching_flag = 1;
				}
				if(fbuf->done_multiple_flag)
					fbuf->done_matching_flag = 1;
			}
		
	}
	//TODO reader may start next dtf transfer session immediately and it may mix up with 
	//with previous transfer session is writer is slow to receive notification that 
	//previous transfer has finished. Probably have to put barrier and make reader wait for confirmation? 
	//but this is problem only if repeadet read and write by same components....
done_match:
	/*Reset flags*/
    if(fbuf->my_mst_info->my_mst == gl_proc.myrank)
		fbuf->my_mst_info->nread_completed = 0;
	fbuf->done_matching_flag = 0;
	fbuf->is_transferring = 0;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished match ioreqs for %s", fbuf->file_path);

	//if(gl_proc.conf.sclltkf_flag)
	DTF_DBG(VERBOSE_ERROR_LEVEL, "dtf_time transfer for %s: %.4f", fbuf->file_path, MPI_Wtime() - t_start);
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
	
	if(gl_proc.msgbuf == NULL) gl_proc.msgbuf = dtf_malloc(gl_proc.conf.data_msg_size_limit);
                    
    sbuf = (unsigned char*)gl_proc.msgbuf; //dtf_malloc(gl_proc.conf.data_msg_size_limit);//
    sbufsz = (int)gl_proc.conf.data_msg_size_limit;

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
					too small (current value %d). Please increase by setting up DFT_DATA_MSG_SIZE_LIMIT", gl_proc.conf.data_msg_size_limit);
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
        if(gl_proc.conf.do_checksum){
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
				int flag;
                assert(sofft != MAX_FILE_NAME );
                DTF_DBG(VERBOSE_DBG_LEVEL, "Send msg to %d, size %d", rdr_rank,(int)sofft);
                /*Send this message*/
                MPI_Request req;
                int dsz = (int)sofft;
                double t_start_comm = MPI_Wtime();
                err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_proc.comps[fbuf->reader_id].comm, &req);
                CHECK_MPI(err);
                
				while(1){
					err = MPI_Iprobe(MPI_ANY_SOURCE, COMP_FINALIZED_TAG, gl_proc.comps[fbuf->reader_id].comm, &flag, MPI_STATUS_IGNORE);
					CHECK_MPI(err);
					
					if(flag){
						err = MPI_Cancel(&req);
						CHECK_MPI(err);
						goto fn_exit;
					}
					err = MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
					CHECK_MPI(err);
					if(flag) break;
				}

                gl_proc.stats_info.t_comm += MPI_Wtime() - t_start_comm;
                gl_proc.stats_info.ndata_msg_sent++;
                gl_proc.stats_info.data_msg_sz += sofft;

                sofft = 0;
            } else {                
                int type_mismatch = 0;
                int cnt_mismatch = 0;
                
                DTF_DBG(VERBOSE_DBG_LEVEL, "Will copy subblock (strt->cnt):");
                for(i=0; i < var->ndims; i++)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "%lld\t --> %lld \t orig: \t %lld\t --> %lld)", new_start[i], new_count[i], start[i], count[i]);              
                
               /*copy the subblock to the sbuf*/
                /*save var_id, start[], count[]*/
                *(MPI_Offset*)(sbuf+sofft) = (MPI_Offset)var_id;
                sofft += sizeof(MPI_Offset);
                memcpy(sbuf+sofft, new_start, var->ndims*sizeof(MPI_Offset));
                sofft += var->ndims*sizeof(MPI_Offset);
                memcpy(sbuf+sofft, new_count, var->ndims*sizeof(MPI_Offset));
                sofft += var->ndims*sizeof(MPI_Offset);
                
				if(var->dtype != ioreq->dtype) type_mismatch = 1;

                DTF_DBG(VERBOSE_DBG_LEVEL, "total data to getput %lld (of nelems %lld)", fit_nelems, nelems);
                
                /*Copy data*/
				for(i = 0; i < var->ndims; i++) 
					if( (new_count[i] != 1) && (new_count[i] < ioreq->count[i])){
						cnt_mismatch = 1;
						break;
					}
					
				if(ioreq->derived_params != NULL) cnt_mismatch = 1;

				if(cnt_mismatch == 0){//it's a continious block of memory
					int req_el_sz;
					MPI_Offset start_cpy_offt;
					MPI_Type_size(ioreq->dtype, &req_el_sz);
					start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, new_start)*req_el_sz;
					DTF_DBG(VERBOSE_DBG_LEVEL, "cpy offt %lld", start_cpy_offt);
					if(type_mismatch){
						convertcpy(ioreq->dtype, var->dtype, (unsigned char*)ioreq->user_buf+start_cpy_offt,
								   (void*)(sbuf+sofft), (int)fit_nelems);
					}else{
						memcpy(sbuf+sofft, (unsigned char*)ioreq->user_buf+start_cpy_offt,fit_nelems*def_el_sz);
					}
					
					if(gl_proc.conf.do_checksum){
						double chsum = compute_checksum(sbuf+sofft, var->ndims, new_count, var->dtype);
						DTF_DBG(VERBOSE_DBG_LEVEL, "chsum contin for req %lu is %.4f", ioreq->id, chsum);
					}
				} else {
					//adjust start coordinate with respect to the user buffer
					MPI_Offset *rel_start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
					for(i = 0; i < var->ndims; i++)
						rel_start[i] = new_start[i] - ioreq->start[i];
					
					if(ioreq->derived_params == NULL)
						get_put_data(var->ndims, ioreq->dtype, var->dtype, ioreq->count, ioreq->user_buf, 
									 rel_start, new_count, sbuf+sofft, DTF_READ, type_mismatch);
					else{
						DTF_DBG(VERBOSE_DBG_LEVEL, "Will extract data from derived data type");
						//adjust with respect to shift of subarray inside to original array
						for(i = 0; i < var->ndims; i++)
							rel_start[i] = ioreq->derived_params->orig_start[i] + rel_start[i];
							
						get_put_data(var->ndims, ioreq->dtype, var->dtype, 
									ioreq->derived_params->orig_array_size, ioreq->user_buf, rel_start, new_count, sbuf+sofft, DTF_READ, type_mismatch);
					}
					dtf_free(rel_start, var->ndims*sizeof(MPI_Offset));
					
					if(gl_proc.conf.do_checksum){
						double chsum = compute_checksum(sbuf+sofft, var->ndims, new_count, var->dtype);
						DTF_DBG(VERBOSE_DBG_LEVEL, "chsum getput for req %lu is %.4f", ioreq->id, chsum);
					}
				}
                
				nblocks_written++;
                gl_proc.stats_info.nbl++;
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
		int flag;
        /*Send the last message*/
        DTF_DBG(VERBOSE_DBG_LEVEL, "Send (last) msg to %d size %d", rdr_rank,(int)sofft);
        MPI_Request req;
  
        int dsz = (int)sofft;
        double t_start_comm = MPI_Wtime();
		err = MPI_Isend((void*)sbuf, dsz, MPI_BYTE, rdr_rank, IO_DATA_TAG, gl_proc.comps[fbuf->reader_id].comm, &req);
		CHECK_MPI(err); 
		while(1){
			err = MPI_Iprobe(MPI_ANY_SOURCE, COMP_FINALIZED_TAG, gl_proc.comps[fbuf->reader_id].comm, &flag, MPI_STATUS_IGNORE);
			CHECK_MPI(err);
			
			if(flag){
				err = MPI_Cancel(&req);
				CHECK_MPI(err);
				goto fn_exit;
			}
			err = MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
			CHECK_MPI(err);
			if(flag) break;
		}
        gl_proc.stats_info.t_comm += MPI_Wtime() - t_start_comm;
        gl_proc.stats_info.ndata_msg_sent++;
        gl_proc.stats_info.data_msg_sz += sofft;
    }

	/*Record the pattern if needed*/
	record_io_pat(fbuf->file_path, rdr_rank, buf, bufsz, fbuf->cur_transfer_epoch);
fn_exit:
    free(new_start);
    free(new_count);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished sending the data");
   // dtf_free(sbuf, sbufsz); //comment out since use gl_proc.msgbuf
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
	int cnt_mismatch;
	
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
   
        /*Find the ioreq that has this block*/
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
					for(i = 0; i < var->ndims; i++)
						assert(start[i]+count[i]<=ioreq->start[i]+ioreq->count[i]);
                    break;
                }
            }
            ioreq = ioreq->next;
        }
        assert(ioreq != NULL);
		
		type_mismatch = 0;
        /*Copy data*/
		MPI_Type_size(var->dtype, &def_el_sz);
		MPI_Type_size(ioreq->dtype, &req_el_sz);
		if(ioreq->dtype != var->dtype) type_mismatch = 1;

        cnt_mismatch = 0;
        for(i = 0; i < var->ndims; i++) 
			if( (count[i] != 1) && (count[i] < ioreq->count[i])){
				cnt_mismatch = 1;
				break;
			}   
        if(ioreq->derived_params != NULL) cnt_mismatch = 1;

        /*Copy data*/
        nelems = 1;
        DTF_DBG(VERBOSE_DBG_LEVEL, "Getput data:");
        for(i = 0; i < var->ndims; i++){
            nelems *= count[i];
        }
        assert(nelems <= (int)(ioreq->req_data_sz - ioreq->get_sz)/req_el_sz);

        DTF_DBG(VERBOSE_DBG_LEVEL, "Will get %d elems for var %d", nelems, var->id);
        assert(ioreq->user_buf != NULL);
        if(cnt_mismatch == 0 ){ //continious block of memory
            MPI_Offset start_cpy_offt = to_1d_index(var->ndims, ioreq->start, ioreq->count, start) * req_el_sz;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Start cpy offt %d", (int)start_cpy_offt);
            
            if(type_mismatch){
                convertcpy(var->dtype, ioreq->dtype, data, (unsigned char*)ioreq->user_buf + start_cpy_offt, nelems);
            } else{
                memcpy((unsigned char*)ioreq->user_buf + start_cpy_offt, data, nelems*def_el_sz);
			}
        } else{
			MPI_Offset *rel_start = dtf_malloc(var->ndims*sizeof(MPI_Offset));
			for(i = 0; i < var->ndims; i++)
				rel_start[i] = start[i] - ioreq->start[i];
			
			if(ioreq->derived_params == NULL)
				get_put_data(var->ndims, var->dtype, ioreq->dtype, ioreq->count, ioreq->user_buf, 
							 rel_start, count, data, DTF_WRITE, type_mismatch);
			else{
				DTF_DBG(VERBOSE_DBG_LEVEL, "Will extract data to derived data type");
				//adjust with respect to shift of subarray inside to original array
				for(i = 0; i < var->ndims; i++)
					rel_start[i] = ioreq->derived_params->orig_start[i] + rel_start[i];
					
				get_put_data(var->ndims, var->dtype, ioreq->dtype,  
							ioreq->derived_params->orig_array_size, ioreq->user_buf, rel_start, 
							count, data, DTF_WRITE, type_mismatch);
			}
			dtf_free(rel_start, var->ndims*sizeof(MPI_Offset));
        }
        gl_proc.stats_info.nbl++;
        offt += nelems*def_el_sz +(nelems*def_el_sz)%sizeof(MPI_Offset);
        ioreq->get_sz += (MPI_Offset)(nelems*req_el_sz);
        
        DTF_DBG(VERBOSE_DBG_LEVEL, "req %lu, var %d, Got %d (expect %d)", ioreq->id, var_id, (int)ioreq->get_sz, (int)ioreq->req_data_sz);
        assert(ioreq->get_sz<=ioreq->req_data_sz);
        
        if(ioreq->get_sz == ioreq->req_data_sz){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Complete req %lu (left %lu), var %d, ", ioreq->id, fbuf->rreq_cnt-1, var->id);
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

            if(gl_proc.conf.do_checksum){
                double chsum = compute_checksum(ioreq->user_buf, var->ndims, ioreq->count, ioreq->dtype);
                DTF_DBG(VERBOSE_DBG_LEVEL, "chsum for req %lu is %.4f", ioreq->id, chsum);
            }
            //delete this ioreq
            delete_ioreq(fbuf, var_id, &ioreq);

            if(fbuf->rreq_cnt == 0){
				fname_pattern_t *pat = find_fname_pattern(fbuf->file_path);
				assert(pat != NULL);
				if(pat->replay_io && pat->rdr_recorded == IO_PATTERN_RECORDED){
					DTF_DBG(VERBOSE_DBG_LEVEL, "Completed rreqs for file %s, won't notify writer because replaying I/O", fbuf->file_path);
					fbuf->done_matching_flag = 1;
				} else {
		
					int err;
					DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE:Completed all rreqs for file %s", fbuf->file_path);
					dtf_msg_t *msg = new_dtf_msg(NULL, 0, DTF_UNDEFINED, READ_DONE_TAG, 1);
					err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->my_mst_info->my_mst,
							  READ_DONE_TAG, gl_proc.comps[gl_proc.my_comp].comm, msg->reqs);
					CHECK_MPI(err);

					//Reader ranks except for root complete the send immediately
					if(gl_proc.myrank != fbuf->my_mst_info->my_mst){
						err = MPI_Wait(msg->reqs, MPI_STATUS_IGNORE);
						CHECK_MPI(err);
						dtf_free(msg, sizeof(dtf_msg_t));
						fbuf->done_matching_flag = 1;
					} else
						ENQUEUE_ITEM(msg, gl_proc.comps[gl_proc.my_comp].out_msg_q);
				}
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
    dtf_msg_t *msg = new_dtf_msg(sbuf, sbuf_sz, DTF_UNDEFINED, FILE_INFO_TAG, 1);
    err = MPI_Isend(sbuf, (int)sbuf_sz, MPI_BYTE, reader_root, FILE_INFO_TAG, gl_proc.comps[fbuf->reader_id].comm, msg->reqs);
    CHECK_MPI(err);
    err = MPI_Test(msg->reqs, &flag, MPI_STATUS_IGNORE);
	CHECK_MPI(err);
	ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->reader_id].out_msg_q);
    gl_proc.stats_info.t_comm += MPI_Wtime() - t_start;
}


static void notify_processes(int group, file_buffer_t *fbuf, void *buf, int bufsz, int msgtag)
{
	int i, err;    
	dtf_msg_t *msg;
	int nprocs; 
	MPI_Datatype dtype;
	
	nprocs = (group == DTF_GROUP_MST) ? fbuf->my_mst_info->nmasters :fbuf->my_mst_info->my_wg_sz;	
	dtype = (bufsz == 0) ? MPI_CHAR : MPI_BYTE;
	
	
	if(bufsz == 0){
		char *fnm = dtf_malloc(MAX_FILE_NAME);
		strcpy(fnm, fbuf->file_path);
		
		msg = new_dtf_msg(fnm, MAX_FILE_NAME, DTF_UNDEFINED, msgtag, nprocs);
	} else
		msg = new_dtf_msg(buf, bufsz, DTF_UNDEFINED, msgtag, nprocs);
	
	msg->reqs[0] = MPI_REQUEST_NULL; //first message is always to itself
	
	if(group == DTF_GROUP_MST)
		for(i = 1; i < nprocs; i++) {	
			err = MPI_Isend(msg->buf, msg->bufsz, dtype, fbuf->my_mst_info->masters[i], msgtag, 
							gl_proc.comps[gl_proc.my_comp].comm, &(msg->reqs[i]));
			CHECK_MPI(err);
		}
	else{
		int rank;
		MPI_Comm_rank(fbuf->comm, &rank);
		int *ranks_file = dtf_malloc(nprocs*sizeof(int));  
		int *ranks_glob = dtf_malloc(nprocs*sizeof(int));

		for(i = 0; i < nprocs; i++)
			ranks_file[i] = rank + i;
			
		translate_ranks(ranks_file, nprocs, fbuf->comm, gl_proc.comps[gl_proc.my_comp].comm, ranks_glob);
		dtf_free(ranks_file, nprocs*sizeof(int));
		
		for(i = 1; i < nprocs; i++) {	
			err = MPI_Isend(msg->buf, msg->bufsz, dtype, ranks_glob[i], msgtag, 
							gl_proc.comps[gl_proc.my_comp].comm, &(msg->reqs[i]));
			CHECK_MPI(err);
		}
		dtf_free(ranks_glob, nprocs*sizeof(int));
	}
	err = MPI_Waitall(nprocs, msg->reqs, MPI_STATUSES_IGNORE);
	CHECK_MPI(err);
	delete_dtf_msg(msg);
    //ENQUEUE_ITEM(msg, gl_proc.comps[gl_proc.my_comp].out_msg_q);
	//progress_send_queue();
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
        case FILE_INFO_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_TAG from %d", src);
            break;
        case FILE_INFO_REQ_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag FILE_INFO_REQ_TAG from %d", src);
            break;
        case MATCH_DONE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag MATCH_DONE_TAG from %d", src);
            break;
        case DONE_MULTIPLE_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag DONE_MULTIPLE_TAG from %d", src);
            break;
        case IO_REQS_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_REQS_TAG from %d", src);
            break;
        case IO_DATA_TAG:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag IO_DATA_TAG from %d", src);
            break;
		case COMP_FINALIZED_TAG:
			DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag COMP_FINALIZED_TAG from %d", src);
			break;	
        default:
            DTF_DBG(VERBOSE_DBG_LEVEL, "Received tag unknown %d from %d", tag, src);
            assert(0);
    }
}

int parse_msg(int comp, int src, int tag, void *rbuf, int bufsz, int is_queued)
{
	MPI_Status status;
	file_buffer_t *fbuf;
	size_t offt = 0;
	char filename[MAX_FILE_NAME];
	int ret = 1;
	int free_buf = 1;
	fname_pattern_t *pat;
	
	if(tag == COMP_FINALIZED_TAG){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Comp %s finalized", gl_proc.comps[comp].name);
		gl_proc.comps[comp].finalized = 1;
		return 0;
	}
	memcpy(filename, rbuf, MAX_FILE_NAME);
	offt += MAX_FILE_NAME ;
	fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
	if(fbuf == NULL) goto fn_exit;
	if(fbuf->iomode==DTF_IO_MODE_MEMORY && !fbuf->is_defined) goto fn_exit;
	
	switch(tag){
			case FILE_INFO_REQ_TAG:
				if(!gl_proc.conf.sclltkf_flag){
					//parse mst info of the other component
					assert(fbuf->cpl_mst_info->nmasters == 0);
					memcpy(&(fbuf->cpl_mst_info->comm_sz), (unsigned char*)rbuf+offt, sizeof(int));
					offt+=sizeof(int);
					memcpy(&(fbuf->cpl_mst_info->nmasters), (unsigned char*)rbuf+offt, sizeof(int));
					offt+=sizeof(int);
					assert(fbuf->cpl_mst_info->nmasters > 0);
					fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
					memcpy(fbuf->cpl_mst_info->masters, (unsigned char*)rbuf+offt, fbuf->cpl_mst_info->nmasters*sizeof(int));
					fbuf->root_reader = fbuf->cpl_mst_info->masters[0];
				}
				
				if(fbuf->my_mst_info->my_mst == gl_proc.myrank){
					
					if(fbuf->root_writer == gl_proc.myrank && fbuf->iomode == DTF_IO_MODE_MEMORY){
						DTF_DBG(VERBOSE_DBG_LEVEL, "I am root writer for file %s, process the file info request, root reader %d",
							fbuf->file_path, fbuf->cpl_mst_info->masters[0] );
						send_file_info(fbuf, fbuf->root_reader);
						//forward info to others
						void *tmp = dtf_malloc(bufsz);
						memcpy(tmp, rbuf, bufsz);
						notify_processes(DTF_GROUP_MST, fbuf, tmp, bufsz, FILE_INFO_REQ_TAG);
					}

					void *tmp2 = dtf_malloc(bufsz);
					memcpy(tmp2, rbuf, bufsz);
					notify_processes(DTF_GROUP_WG, fbuf, tmp2, bufsz, FILE_INFO_REQ_TAG);					
					
				} else 
					DTF_DBG(VERBOSE_DBG_LEVEL, "Got info about the other component from master");
				
				pat = find_fname_pattern(fbuf->file_path);
				assert(pat != NULL);
				if(pat->replay_io){
					assert(pat->finfo_sz == 0);
					pat->finfo_sz = bufsz;
					pat->finfo = rbuf;
					free_buf = 0;
				}
				fbuf->cpl_info_shared = 1;			
				break;
			case IO_REQS_TAG:
				if( (comp != gl_proc.my_comp) && (fbuf->cpl_mst_info->comm_sz == 0)) goto fn_exit;
				if( gl_proc.comps[comp].finalized){
						DTF_DBG(VERBOSE_DBG_LEVEL, "Discard message as the component has started finalizing.");
						goto discard;
				}
				if(!fbuf->is_transferring) goto fn_exit;  //Not ready to process this message yet
				DTF_DBG(VERBOSE_DBG_LEVEL, "Received reqs from %d, comp %d (my comp %d), bufsz %d", src, comp,gl_proc.my_comp, (int)bufsz);
				parse_ioreqs(fbuf, (unsigned char*)rbuf+offt, bufsz - offt, src, gl_proc.comps[comp].comm);
				break;
			case IO_DATA_REQ_TAG:
				if(!fbuf->is_transferring) goto fn_exit;  //Not ready to process this message yet
				send_data(fbuf, (unsigned char*)rbuf+offt, bufsz - offt);
				break;
			case IO_DATA_TAG:
				if(!fbuf->is_transferring) goto fn_exit;  //Not ready to process this message yet
				recv_data_rdr(fbuf, (unsigned char*)rbuf + offt, bufsz - offt);
				break;
			case READ_DONE_TAG:
				fbuf->my_mst_info->nread_completed++;
				DTF_DBG(VERBOSE_DBG_LEVEL, "Recv read done for file %s from %d in comp %d (tot %d)", fbuf->file_path,
						src, comp, fbuf->my_mst_info->nread_completed);

				if(gl_proc.my_comp == fbuf->writer_id){
					assert(comp == fbuf->reader_id);
					//assert(gl_proc.myrank == fbuf->root_writer);
					if(fbuf->my_mst_info->nread_completed == fbuf->cpl_mst_info->nmasters){
						//notify_processes(DTF_GROUP_MST, fbuf, NULL, 0, MATCH_DONE_TAG);
						notify_processes(DTF_GROUP_WG, fbuf,  NULL, 0, MATCH_DONE_TAG);
						fbuf->done_matching_flag = 1;
					}
					
					
				} else { /*reader comp*/
					
					assert(fbuf->my_mst_info->my_mst == gl_proc.myrank);
					
					if(fbuf->my_mst_info->nread_completed == fbuf->my_mst_info->my_wg_sz){
							int err, i;
							char *fname = dtf_malloc(MAX_FILE_NAME);
							strcpy(fname, fbuf->file_path);
							DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer masters %d that my workgroup completed reading", fbuf->root_writer);
							
							dtf_msg_t *msg = new_dtf_msg(fname, MAX_FILE_NAME, DTF_UNDEFINED, READ_DONE_TAG, fbuf->cpl_mst_info->nmasters);
							for(i=0; i < fbuf->cpl_mst_info->nmasters; i++){
								err = MPI_Isend(fname, MAX_FILE_NAME, MPI_CHAR, fbuf->cpl_mst_info->masters[i], READ_DONE_TAG, 
												gl_proc.comps[fbuf->writer_id].comm, &(msg->reqs[i]));
								CHECK_MPI(err);
							}
							//~ err = MPI_Test(msg->reqs, &flag, MPI_STATUS_IGNORE);
							err = MPI_Waitall(fbuf->cpl_mst_info->nmasters, msg->reqs, MPI_STATUSES_IGNORE);
							CHECK_MPI(err);
							delete_dtf_msg(msg);
							//~ ENQUEUE_ITEM(msg, gl_proc.comps[fbuf->writer_id].out_msg_q);
							fbuf->done_matching_flag = 1;
							DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
					}
				}
				break;
			case MATCH_DONE_TAG:
				if(fbuf->my_mst_info->my_mst == gl_proc.myrank)
					notify_processes(DTF_GROUP_WG, fbuf, NULL, 0, MATCH_DONE_TAG);
				DTF_DBG(VERBOSE_DBG_LEVEL, "Done matching flag set for file %s", fbuf->file_path);
				fbuf->done_matching_flag = 1;
				break;
			case DONE_MULTIPLE_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL, "Recv done multiple tag for %s", fbuf->file_path);

				if(fbuf->my_mst_info->my_mst == gl_proc.myrank){
					DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writers that multiple matching for %s completed", fbuf->file_path);
					/*Tell other writer ranks that they can complete matching*/
					notify_processes(DTF_GROUP_WG, fbuf, NULL, 0, DONE_MULTIPLE_TAG);
				}
				fbuf->done_multiple_flag = 1;
				break;
		   case FILE_READY_TAG:
				DTF_DBG(VERBOSE_DBG_LEVEL,   "Process FILE_READY notif for %s", (char*)rbuf);
				fbuf->is_ready = 1;
				break;
		   default:
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: unknown tag %d", status.MPI_TAG);
				assert(0);
		}

discard:
	if(free_buf)dtf_free(rbuf, bufsz);	
	return ret;
	
fn_exit:
	if(!is_queued){
		DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: can't process message for file %s yet, queue it to process later", filename); 
		dtf_msg_t *msg = new_dtf_msg(rbuf, bufsz, src, tag, 0);
		ENQUEUE_ITEM(msg, gl_proc.comps[comp].in_msg_q);
	}
	return 0;
}

void progress_comm()
{
    MPI_Status status;
    int comp, flag, err, bufsz, tag, src;
    gl_proc.stats_info.nprogress_call++;
    double t_st, t_start_comm;
    void *rbuf = NULL;
	t_st = MPI_Wtime();
	progress_recv_queue();
    progress_send_queue();

    for(comp = 0; comp < gl_proc.ncomps; comp++){
        if( (gl_proc.comps[comp].comm == MPI_COMM_NULL) || (comp == gl_proc.my_comp)){
            continue;
        }

        while(1){
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_proc.comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                gl_proc.stats_info.idle_time += MPI_Wtime() - t_st;
                if(MPI_Wtime() - gl_proc.stats_info.t_idle > DTF_TIMEOUT){
					DTF_DBG(VERBOSE_ERROR_LEVEL, "Process have been idle for %.2f seconds. Consider that it hang and abort execution", MPI_Wtime() - gl_proc.stats_info.t_idle);
					MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
				}
                break;
            } 
            gl_proc.stats_info.t_idle = MPI_Wtime();
            tag = status.MPI_TAG;
            src = status.MPI_SOURCE;
			MPI_Get_count(&status, MPI_BYTE, &bufsz);
			if(bufsz > 0) rbuf = dtf_malloc(bufsz);
            print_recv_msg(tag, src);
            t_start_comm = MPI_Wtime();
			err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, tag, gl_proc.comps[comp].comm, &status);
			gl_proc.stats_info.t_comm += MPI_Wtime() - t_start_comm;
			CHECK_MPI(err);
			parse_msg(comp, src, tag, rbuf, bufsz, 0);
        }
    }
    rbuf = NULL;
    //check messages in my component
    comp = gl_proc.my_comp;
    while(1){
            err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_proc.comps[comp].comm, &flag, &status);
            CHECK_MPI(err);

            if(!flag){
                gl_proc.stats_info.idle_time += MPI_Wtime() - t_st;
                break;
            }         
			tag = status.MPI_TAG;
			src = status.MPI_SOURCE;
			MPI_Get_count(&status, MPI_BYTE, &bufsz);
			if(bufsz > 0) rbuf = dtf_malloc(bufsz);
            print_recv_msg(tag, src);
            t_start_comm = MPI_Wtime();
			err = MPI_Recv(rbuf, bufsz, MPI_BYTE, src, tag, gl_proc.comps[comp].comm, &status);
			gl_proc.stats_info.t_comm += MPI_Wtime() - t_start_comm;
			CHECK_MPI(err);
			parse_msg(comp, src, tag, rbuf, bufsz, 0);
	}
    gl_proc.stats_info.t_progr_comm += MPI_Wtime() - t_st;
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
	req->var_id = varid;
	req->ndims = ndims;
	req->next = NULL;
	req->rw_flag = rw_flag;
	req->dtype = dtype;

    MPI_Type_size(dtype, &el_sz);

	if(ndims > 0){
        int i;
        req->req_data_sz = el_sz;
        for(i=0;i<ndims;i++)
            req->req_data_sz *= count[i];
        
        req->start = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        memcpy((void*)req->start, (void*)start, sizeof(MPI_Offset)*ndims);
        req->count = (MPI_Offset*)dtf_malloc(sizeof(MPI_Offset)*ndims);
        memcpy((void*)req->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else{
        req->start = NULL;
        req->count = NULL;
        req->req_data_sz = el_sz;
    }
    req->user_buf = buf;
    /*checksum*/
    //~ if( (rw_flag == DTF_WRITE) && gl_proc.conf.do_checksum && (dtype == MPI_DOUBLE || dtype == MPI_FLOAT)){
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

void notify_complete_multiple(file_buffer_t *fbuf)
{
    int i, err;
    
    if(fbuf->root_reader == gl_proc.myrank){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Will notify writer masters that completed multiple for %s", fbuf->file_path);
		char *fnm = dtf_malloc(MAX_FILE_NAME);
		strcpy(fnm, fbuf->file_path);
		dtf_msg_t *msg = new_dtf_msg(fnm, MAX_FILE_NAME, DTF_UNDEFINED, DONE_MULTIPLE_TAG, fbuf->cpl_mst_info->nmasters);
            
        for(i = 0; i < fbuf->cpl_mst_info->nmasters; i++){
            err = MPI_Isend(fnm,MAX_FILE_NAME, MPI_CHAR, fbuf->cpl_mst_info->masters[i], DONE_MULTIPLE_TAG, gl_proc.comps[fbuf->writer_id].comm, &(msg->reqs[i]));
            CHECK_MPI(err);
            
        }
        ENQUEUE_ITEM(msg, gl_proc.comps[gl_proc.my_comp].out_msg_q);
		progress_send_queue();
    }
}
