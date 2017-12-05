#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "dtf_util.h"
#include "dtf_req_match.h"
#include "dtf_file_buffer.h"
#include "dtf.h"

/*API for handling rb_tree in write_db_item*/

/*
int var_cmp(const void *a, const void *b)
{
  if( ((dtf_var_t*)a)->id > ((dtf_var_t*)b)->id ) return 1;
  if( ((dtf_var_t*)a)->id < ((dtf_var_t*)b)->id ) return -1;
  return 0;
}

void var_print(const void *var)
{
  DTF_DBG(VERBOSE_DBG_LEVEL, "(varid %d)", ((dtf_var_t*)var)->id);
} */



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

        if(ptr->slink_name!=NULL &&  (strstr(file_path, ptr->slink_name)!=NULL || strstr(ptr->slink_name, file_path)!=NULL) )
           break;

        ptr = ptr->next;
    }
		
	//DTF_DBG(VERBOSE_DBG_LEVEL, "Return fbuf %p", (void*)ptr);
    return ptr;
}

void delete_var(dtf_var_t* var)
{
	assert(var->ioreqs == NULL);
    dtf_free(var->shape, var->ndims*sizeof(MPI_Offset));
    dtf_free(var, sizeof(dtf_var_t));
}

void add_var(file_buffer_t *fbuf, dtf_var_t *var)
{
    //rb_red_blk_node *var_node = RBTreeInsert(fbuf->vars, var, 0);
    //assert(var_node != NULL);
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

    if(fbuf == NULL)
        return;

    assert(gl_filebuf_list != NULL);

    if(fbuf->header != NULL)
        dtf_free(fbuf->header, fbuf->hdr_sz);

	{
		int is_scale = 0;
		
		char *c = getenv("DTF_SCALE");
		if(c != NULL)
			is_scale = atoi(c);
			
		if(is_scale){
			if(fbuf->mst_info->is_master_flag)
				clean_iodb(fbuf->mst_info->iodb, fbuf->nvars);
		
			delete_ioreqs(fbuf,1);
		}
	}

    //RBTreeDestroy(fbuf->vars);
    for(i = 0; i < fbuf->nvars; i++)
        delete_var(fbuf->vars[i]);
    dtf_free(fbuf->vars, fbuf->nvars*sizeof(dtf_var_t*));
	
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
    
    if(fbuf->mst_info->iodb != NULL){
        finalize_iodb(fbuf);
    }
    if(fbuf->mst_info != NULL){
        if(fbuf->mst_info->masters != NULL)
            dtf_free(fbuf->mst_info->masters, fbuf->mst_info->nmasters*sizeof(int));
        if(fbuf->mst_info->my_wg != NULL)
            dtf_free(fbuf->mst_info->my_wg,(fbuf->mst_info->my_wg_sz - 1)*sizeof(int));
        dtf_free(fbuf->mst_info, sizeof(master_info_t));
    }
    
    if( (fbuf->reader_id == gl_my_comp_id) && (fbuf->root_reader == gl_my_rank) && (fbuf->slink_name != NULL)){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Remove symbolic link %s", fbuf->slink_name);
		int err;
		char *dir = NULL;
		char wdir[MAX_FILE_NAME]="\0";
		char slink[MAX_FILE_NAME]="\0";

		dir = getcwd(wdir, MAX_FILE_NAME);
		assert(dir != NULL);
		sprintf(slink, "%s/%s", wdir, fbuf->slink_name);

		DTF_DBG(VERBOSE_DBG_LEVEL, "Remove symlink %s", slink);
		err = unlink(slink);
		if(err)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: error deleting symbolic link %s", slink);
	}

    
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
    pat->rdr = -1;
    pat->wrt = -1;
    pat->ignore_io = 0;
    pat->slink_name = NULL;
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
	DTF_DBG(VERBOSE_DBG_LEVEL, "Matched against pattern %s", pat->fname);
	buf = dtf_malloc(sizeof(struct file_buffer));
    assert( buf != NULL );
    buf->next = NULL;
    buf->prev = NULL;
    buf->is_ready = 0;   //TODO what to do with this flag when letkf openes the file the second time?
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
    buf->done_match_confirm_flag = DTF_UNDEFINED;
    buf->comm = MPI_COMM_NULL;
    buf->root_writer = -1;
    buf->root_reader = -1;
    
    buf->nwriters = 0;
    buf->mst_info = dtf_malloc(sizeof(master_info_t));
    assert(buf->mst_info != NULL);
    buf->mst_info->masters = NULL;
    buf->mst_info->is_master_flag = 0;
    buf->mst_info->my_wg_sz = 0;
    buf->mst_info->my_wg = NULL;
    buf->mst_info->nmasters = 0;
    buf->mst_info->iodb = NULL;
    buf->mst_info->nrranks_completed = 0;
    buf->mst_info->nranks_opened = 0;
    buf->is_matching_flag = 0;
	strcpy(buf->file_path, file_path);
	buf->reader_id = pat->rdr;
	buf->writer_id = pat->wrt;
	if(pat->slink_name != NULL){
		buf->slink_name = dtf_malloc(MAX_FILE_NAME*sizeof(char));
		assert(buf->slink_name != NULL);
		strcpy(buf->slink_name, pat->slink_name);
	} else 
		buf->slink_name = NULL;
	buf->iomode = pat->iomode;
	buf->ignore_io = pat->ignore_io;
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
        for(i=0; i<ndims;i++)
            var->shape[i] = shape[i];
        var->max_dim = shape[ndims-1] > shape[0] ? ndims-1 : 0;
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


int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start, const MPI_Offset *count )
{
    dtf_var_t *var = fbuf->vars[varid];

    if(var->ndims > 0){
        int i;

      /*  if(frt_indexing){
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[var->ndims-i-1]){
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[var->ndims-i-1]);
                    return 1;
                }
        } else { */
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[i]){
                            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[i]);
                            return 1;
                }
       /* } */
    }
//    else {
//        if( (start != NULL) || (count != NULL)){
//            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d is a scalar variable but trying to read an array", varid);
//            return 1;
//        }
//    }
    return 0;
}
