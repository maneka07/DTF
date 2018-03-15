#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "dtf_util.h"
#include "dtf_req_match.h"
#include "dtf_file_buffer.h"
#include "dtf.h"

static int match_str(char* pattern, const char* filename)
{
    int ret = 0;
    char *next_token, *cur_token;
    int next_pos, cur_pos, token_len;
    char *match_str;
	
    if(strlen(pattern)== 0 || strlen(filename)==0)
        return ret;
        
	next_token = strchr(pattern, '%');
	
	if(next_token == NULL){
		//not a name pattern, just check for inclusion
		if(strstr(pattern, filename)!=NULL || strstr(filename, pattern)!=NULL)
			ret = 1;
	} else {
		cur_pos = 0;
		cur_token = &pattern[cur_pos];
		
		while( next_token != NULL){
			next_pos = next_token-pattern;
			token_len = next_token - cur_token;
			pattern[next_pos] = '\0';
			
			if(next_pos - cur_pos > 0){ //check for inclusion
				match_str = strstr(filename, cur_token);
				ret = (match_str == NULL) ? 0 : 1;
			}
			pattern[next_pos] = '%';
			
			if(!ret) break;
			
			cur_pos = next_pos+1;
			cur_token = &pattern[cur_pos];
			filename = match_str + token_len;
			next_token = strchr(cur_token, '%');
		}
		if( ret && (cur_pos != strlen(pattern)))
			//check the last token
			ret = (strstr(filename, cur_token) == NULL) ? 0 : 1;
	}
	
    return ret;
}

static void init_mst_info(master_info_t* mst_info)
{
    mst_info->masters = NULL;
    mst_info->my_mst = -1;
    mst_info->my_wg_sz = 0;
    mst_info->my_wg = NULL;   
    mst_info->nmasters = 0;
    mst_info->iodb = NULL;
    mst_info->nread_completed = 0;
}

/***************************File structure**************************************/

file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid)
{
	/* If a file buffer for this file exists, return it. Otherwise, check if file name
	 * matches any file name pattern. If yes, create a new file buffer and return it.
	 * Otherwise, return NULL*/
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
       if((ncid >= 0) && (ptr->ncid == ncid))
                break;
       if( (file_path != NULL) &&
           (strlen(file_path)!=0 && ( (strstr(file_path, ptr->file_path)!=NULL || strstr(ptr->file_path, file_path)!=NULL)))){
		   //DTF_DBG(VERBOSE_DBG_LEVEL, "Found by name in list %s", ptr->file_path);
           break;
	   }

       ptr = ptr->next;
    }

	//DTF_DBG(VERBOSE_DBG_LEVEL, "Return fbuf %p", (void*)ptr);
    return ptr;
}

/***************************File structure**************************************/

void delete_var(file_buffer_t *fbuf, dtf_var_t* var)
{
	assert(var->ioreqs == NULL);
    dtf_free(var->shape, var->ndims*sizeof(MPI_Offset));
    dtf_free(var, sizeof(dtf_var_t));
    fbuf->nvars--;
}

void add_var(file_buffer_t *fbuf, dtf_var_t *var)
{
    assert(var->id == fbuf->nvars);
    dtf_var_t **tmp = (dtf_var_t**)realloc(fbuf->vars, (fbuf->nvars+1)*sizeof(dtf_var_t*));
    assert(tmp != NULL);
    fbuf->vars = tmp;
    fbuf->nvars++;
    DTF_DBG(VERBOSE_DBG_LEVEL, "var id %d, cnt %d", var->id, fbuf->nvars);

    fbuf->vars[fbuf->nvars-1] = var;
    gl_stats.malloc_size += sizeof(dtf_var_t*);
}

void delete_file_buffer(file_buffer_t* fbuf)
{
    int i;
    int nvars=fbuf->nvars;

    if(fbuf == NULL)
        return;

    assert(gl_filebuf_list != NULL);

    if(fbuf->header != NULL)
        dtf_free(fbuf->header, fbuf->hdr_sz);

    if(fbuf->my_mst_info != NULL){

        if(fbuf->my_mst_info->iodb != NULL){
            clean_iodb(fbuf->my_mst_info->iodb, fbuf->nvars);
            dtf_free(fbuf->my_mst_info->iodb, sizeof(ioreq_db_t));
        }

        if(fbuf->my_mst_info->masters != NULL)
            dtf_free(fbuf->my_mst_info->masters, fbuf->my_mst_info->nmasters*sizeof(int));
        if(fbuf->my_mst_info->my_wg != NULL)
            dtf_free(fbuf->my_mst_info->my_wg,(fbuf->my_mst_info->my_wg_sz - 1)*sizeof(int));
        dtf_free(fbuf->my_mst_info, sizeof(master_info_t));
    }

    if(fbuf->cpl_mst_info != NULL){
        if(fbuf->cpl_mst_info->masters != NULL)
            dtf_free(fbuf->cpl_mst_info->masters, fbuf->cpl_mst_info->nmasters*sizeof(int));
        assert(fbuf->cpl_mst_info->my_wg == NULL);
        dtf_free(fbuf->cpl_mst_info, sizeof(master_info_t));
    }

    delete_ioreqs(fbuf,1); 

    for(i = 0; i < nvars; i++)
        delete_var(fbuf, fbuf->vars[i]);
    dtf_free(fbuf->vars, nvars*sizeof(dtf_var_t*));

	if(fbuf->ioreq_log != NULL){
		io_req_log_t *ior = fbuf->ioreq_log;
		while(ior != NULL){
			fbuf->ioreq_log = fbuf->ioreq_log->next;

			if(ior->rw_flag == DTF_READ)
				fbuf->rreq_cnt--;
			else
				fbuf->wreq_cnt--;

			if(gl_conf.buffered_req_match && (ior->rw_flag == DTF_WRITE)) dtf_free(ior->user_buf, ior->user_buf_sz);
			if(ior->start != NULL) dtf_free(ior->start, ior->ndims*sizeof(MPI_Offset));
			if(ior->count != NULL)dtf_free(ior->count, ior->ndims*sizeof(MPI_Offset));
			dtf_free(ior, sizeof(io_req_log_t));
			ior = fbuf->ioreq_log;
		}
	}

    assert(fbuf->rreq_cnt == 0);
    assert(fbuf->wreq_cnt == 0);

    if(gl_filebuf_list == fbuf)
		gl_filebuf_list = fbuf->next;

	if(fbuf->prev != NULL)
		fbuf->prev->next = fbuf->next;
	if(fbuf->next != NULL)
		fbuf->next->prev = fbuf->prev;

    dtf_free(fbuf, sizeof(file_buffer_t));
    gl_stats.nfiles--;
}

fname_pattern_t *new_fname_pattern()
{
    fname_pattern_t *pat = dtf_malloc(sizeof(fname_pattern_t));
    assert(pat != NULL);
	pat->fname[0]='\0';
    pat->iomode = DTF_UNDEFINED;
    pat->excl_fnames = NULL;
    pat->nexcls = 0;
    pat->next = NULL;
    pat->comp1 = DTF_UNDEFINED;
    pat->comp2 = DTF_UNDEFINED;
    pat->ignore_io = 0;
    pat->replay_io = DTF_UNDEFINED;
    pat->rdr_recorded = DTF_UNDEFINED;
    pat->wrt_recorded = DTF_UNDEFINED;
    pat->write_only = 0;
    pat->io_pats = NULL;
    return pat;
}



file_buffer_t *create_file_buffer(fname_pattern_t *pat, const char* file_path)
{

	file_buffer_t *buf;

	DTF_DBG(VERBOSE_ALL_LEVEL, "Create file buffer for file %s", file_path);

	if(strlen(file_path) > MAX_FILE_NAME){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: filename %s longer than MAX_FILE_NAME (%d)", file_path, MAX_FILE_NAME);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	buf = dtf_malloc(sizeof(struct file_buffer));
    assert( buf != NULL );
    buf->next = NULL;
    buf->prev = NULL;
    buf->is_ready = 0;
    buf->vars = NULL; //RBTreeCreate(var_cmp, var_destroy, NullFunction, var_print, NullFunction);
    buf->nvars = 0;
    buf->header = NULL;
    buf->hdr_sz = 0;
    buf->rreq_cnt = 0;
    buf->wreq_cnt = 0;
    buf->ioreq_log = NULL;
    buf->ncid = -1;
    buf->done_matching_flag = 0;
    buf->fready_notify_flag = DTF_UNDEFINED;
    buf->sync_comp_flag = 0;
    buf->comm = MPI_COMM_NULL;
    buf->cpl_info_shared = 0;
    buf->my_mst_info = dtf_malloc(sizeof(master_info_t));
    assert(buf->my_mst_info != NULL);
    init_mst_info(buf->my_mst_info);    
	buf->my_mst_info->iodb = dtf_malloc(sizeof(struct ioreq_db));
	assert(buf->my_mst_info->iodb != NULL);
	init_iodb(buf->my_mst_info->iodb);
            
    buf->cpl_mst_info = dtf_malloc(sizeof(master_info_t));
    assert(buf->cpl_mst_info != NULL);
    init_mst_info(buf->cpl_mst_info);
	strcpy(buf->file_path, file_path);
	buf->reader_id = -1;
	buf->writer_id = -1;
	buf->root_writer = -1;
    buf->root_reader = -1;
    //buf->omode = DTF_UNDEFINED;
    buf->is_defined = 0;
	buf->iomode = pat->iomode;
	buf->ignore_io = pat->ignore_io;
	buf->cur_transfer_epoch = 0;
	buf->has_unsent_ioreqs = 0;
	buf->t_last_sent_ioreqs = 0;
	buf->is_transfering = 0;
	//insert
	if(gl_filebuf_list == NULL)
		gl_filebuf_list = buf;
	else{
		buf->next = gl_filebuf_list;
		buf->next->prev = buf;
		gl_filebuf_list = buf;
	}
	gl_stats.nfiles++;
	return buf;
}

dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i;

    dtf_var_t *var = (dtf_var_t*)dtf_malloc(sizeof(dtf_var_t));
    assert(var!=NULL);

    /*Initialize whatever we can initialize at this stage*/
    var->id = varid;
    if(ndims > 0){ //non scalar
        var->shape = (MPI_Offset*)dtf_malloc(ndims*sizeof(MPI_Offset));
        //~ int max_dim = 0;

        for(i=0; i<ndims;i++){
            var->shape[i] = shape[i];
            //~ if(shape[i] > shape[max_dim])
				//~ max_dim = i;
		}

        var->max_dim = 0; //max_dim; //0;
        //~ if(var->shape[0] ==  DTF_UNLIMITED && ndims > 1){
			//~ var->max_dim = 1;
		//~ }
        //~ DTF_DBG(VERBOSE_DBG_LEVEL, "Maxdim for var id %d is %d", varid, var->max_dim);

    } else {
        var->shape = NULL;
		var->max_dim = -1;
	}
    var->ndims = ndims;
    var->dtype = dtype;
    var->checksum = 0;
    var->ioreqs = NULL;
    return var;
}



fname_pattern_t *find_fname_pattern(const char *filename)
{
	int i;
	fname_pattern_t *pat;
	int found = 0;
	
	if (strlen (filename) == 0) return NULL;
		
    pat = gl_fname_ptrns;
	while(pat != NULL){
		
		if(match_str(pat->fname, filename)){
			found = 1;
		 /*First see if the filename matches against any of the exclude patterns*/
			for(i = 0; i < pat->nexcls; i++)
				if(match_str(pat->excl_fnames[i], filename)){
					found = 0;
					break;					
				}
		}
		if(found)
			break;
		pat = pat->next;
	}
	return pat;
}


MPI_Offset read_write_var(file_buffer_t *fbuf,
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
    io_req_t *req;
    int i;
    int def_el_sz, req_el_sz;
    MPI_Offset nelems;
    double t_start = MPI_Wtime();

    if(rw_flag == DTF_READ){
        if(fbuf->reader_id==gl_my_comp_id){
          assert(fbuf->is_ready);
        } else{
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: writer process tries to read file %s (var %d)", fbuf->file_path, varid);
            //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(imap != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    dtf_var_t *var = fbuf->vars[varid];
    DTF_DBG(VERBOSE_DBG_LEVEL, "rw call %d for %s (ncid %d) var %d", rw_flag,fbuf->file_path, fbuf->ncid, var->id);
    for(i = 0; i < var->ndims; i++)
			DTF_DBG(VERBOSE_DBG_LEVEL, "  %lld --> %lld", start[i], count[i]);
    /*check number of elements to read*/
    nelems = 0;
    if(var->ndims == 0)
        nelems = 1;
    else
        if(count != NULL){
            int i;
            nelems = count[0];
            for(i = 1; i < var->ndims; i++)
                nelems *= count[i];

            if(nelems == 0){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Nothing to read or write");
                return 0;
            }
        }
    assert(nelems != 0);

    MPI_Type_size(var->dtype, &def_el_sz);
    MPI_Type_size(dtype, &req_el_sz);

    if(def_el_sz != req_el_sz)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Warning: var %d el_sz mismatch (def %d-bit, access %d).", var->id, def_el_sz, req_el_sz);
    
    //assert(var->dtype == dtype);

	//~ if(rw_flag == DTF_WRITE){
		//~ DTF_DBG(VERBOSE_DBG_LEVEL, "------------WRITE IOREQ--------:");
		
		//~ for(i = 0; i < nelems; i++)
			//~ printf("%.2f\t", ((double*)buf)[i]);
		//~ printf("\n");
	//~ }
        
	/*NOTE: Because dtype may be a derivative MPI type and differ from var->dtype,
	we ignore it. Start and count parameters are supposed to be with respect to
	element size for var->dtype*/
	int buffered = gl_conf.buffered_req_match;

	if(rw_flag == DTF_READ)
		buffered = 0;

	if( gl_scale && (var->ndims <= 1) && (rw_flag == DTF_WRITE))
		 /*This is specifically for SCALE-LETKF since they overwrite the
		  user buffer in every time frame iteration */
		buffered = 1;

	req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, varid, var->ndims, dtype, start, count, buf, rw_flag, buffered);
	
	if( gl_scale && (var->ndims <= 1) && (rw_flag == DTF_WRITE))
		 /*This is specifically for SCALE-LETKF since they overwrite the
		  user buffer in every time frame iteration */
		req->is_permanent = 1; //dont delete this req when cleaning the list of ioreqs

	if(gl_conf.do_checksum && (rw_flag == DTF_WRITE))
		var->checksum += req->checksum;

	if(rw_flag == DTF_READ)
		fbuf->rreq_cnt++;
	else
		fbuf->wreq_cnt++;

	/*Enqueue the request to the head*/
	if(var->ioreqs == NULL)
		var->ioreqs = req;
	else{
		/*Check if some data is overwritten (just to print out a warning message).
		  Becase the new I/O req is pushed to the head of the queue, the
		  writer will access the newest data.*/
		io_req_t *tmpreq = var->ioreqs;
		while(tmpreq != NULL){
			if(req->rw_flag == DTF_WRITE){
				int overlap = 0;
				for(i = 0; i < var->ndims; i++ )
					if( (req->start[i] >= tmpreq->start[i]) && (req->start[i] < tmpreq->start[i] + tmpreq->count[i]))
						overlap++;
					else
						break;

				if(overlap == var->ndims){
					DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: overwriting var %d data: (old (start,count) --> new (start,count)", var->id);
					for(i = 0; i < var->ndims; i++)
						DTF_DBG(VERBOSE_DBG_LEVEL, "(%lld, %lld) --> (%lld, %lld)", tmpreq->start[i], tmpreq->count[i], req->start[i], req->count[i]);
				}
			}
			tmpreq = tmpreq->next;
		}
		var->ioreqs->prev = req;
		req->next = var->ioreqs;
		var->ioreqs = req;
	}
    
    fbuf->has_unsent_ioreqs = 1;
    
    if(MPI_Wtime() - fbuf->t_last_sent_ioreqs >= gl_conf.t_send_ioreqs_freq){
		//Send request to master immediately
		if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
			send_ioreqs_by_var(fbuf);
		else //if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
			send_ioreqs_by_block(fbuf);
	}

    MPI_Type_size(dtype, &el_sz);
    ret = 1;
    for(i = 0; i < var->ndims; i++)
        ret *= count[i];
    ret *= el_sz;

    gl_stats.t_rw_var += MPI_Wtime() - t_start;
    return ret;
}
