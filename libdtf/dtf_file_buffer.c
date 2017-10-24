#include <assert.h>
#include <string.h>
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

//int match_str(char* pattern, char* filename)
//{
//    int ret = 0;
//    if(strlen(pattern)== 0 || strlen(filename)==0)
//        return ret;
//
//    char *substr = strtok(pattern, "...");
//    while(substr != NULL){
//        if(strstr(filename, substr) == NULL){
//            ret = 1;
//        }
//        substr = strtok(NULL, "...");
//    }
//}

file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid)
{
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
       if((ncid >= 0) && (ptr->ncid == ncid))
                break;
       if( (file_path != NULL) &&
           (strlen(file_path)!=0 && ( (strstr(file_path, ptr->file_path)!=NULL || strstr(ptr->file_path, file_path)!=NULL))))
           break;

        if(strlen(ptr->alias_name)!=0 &&  (strstr(file_path, ptr->alias_name)!=NULL || strstr(ptr->alias_name, file_path)!=NULL) )
                break;

        ptr = ptr->next;
    }

    return ptr;
}

void delete_var(dtf_var_t* var)
{
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

void add_file_buffer(file_buffer_t** buflist, file_buffer_t* buf)
{
    if(*buflist == NULL){
        *buflist = buf;
    } else {
        file_buffer_t* tmp = *buflist;
        while(tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = buf;
    }
}

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* fbuf)
{
    int i;

    if(fbuf == NULL)
        return;

    assert(buflist != NULL);

    if(fbuf->header != NULL)
        dtf_free(fbuf->header, fbuf->hdr_sz);

    //RBTreeDestroy(fbuf->vars);
    for(i = 0; i < fbuf->nvars; i++)
        delete_var(fbuf->vars[i]);
    dtf_free(fbuf->vars, fbuf->nvars*sizeof(dtf_var_t*));

    assert(fbuf->ioreqs == NULL);
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
    dtf_free(fbuf, sizeof(file_buffer_t));
}

file_buffer_t* new_file_buffer()
{
    file_buffer_t *buf;

    buf = (file_buffer_t*)dtf_malloc(sizeof(struct file_buffer));
    assert( buf != NULL );
    buf->file_path[0]='\0';
    buf->alias_name[0]='\0';
    buf->slink_name[0]='\0';
    buf->next = NULL;
    buf->reader_id = -1;
    buf->writer_id=-1;
    buf->is_ready = 0;
    buf->vars = NULL; //RBTreeCreate(var_cmp, var_destroy, NullFunction, var_print, NullFunction);
    buf->nvars = 0;
    buf->header = NULL;
    buf->hdr_sz = 0;
    buf->iomode = DTF_UNDEFINED;
    buf->rreq_cnt = 0;
    buf->wreq_cnt = 0;
    buf->ioreqs = NULL;
    buf->ncid = -1;
    buf->explicit_match = 0;
    buf->done_matching_flag = 0;
    buf->rdr_closed_flag = 0;
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

        //if(shape[0] == DTF_UNLIMITED)
          //  var->max_dim = 0;
        //else {
            int maxdim = 0;
            MPI_Offset maxdimlen = shape[0];
            for(i = 1; i < ndims; i++)
                if(shape[i] > maxdimlen){
                    maxdim = i;
                    maxdimlen = shape[i];
                }
            var->max_dim = maxdim;
        //}
    }
    else
        var->shape = NULL;
    var->ndims = ndims;
    var->dtype = dtype;
    var->checksum = 0;
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
