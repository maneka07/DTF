#include <assert.h>
#include <string.h>
#include "dtf_util.h"
#include "dtf_req_match.h"
#include "dtf_file_buffer.h"
#include "dtf.h"
#include "dtf_common.h"


file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid)
{
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
       if((ncid != -1) && (ptr->ncid == ncid))
                break;
       if( (file_path != NULL) &&
           ( (strlen(file_path)!=0 && strstr(file_path, ptr->file_path)!=NULL) ||
             (strlen(ptr->alias_name)!=0 && strstr(file_path, ptr->alias_name)!=NULL) ) )
                break;

        ptr = ptr->next;
    }

    return ptr;
}

dtf_var_t* find_var(dtf_var_t* varlist, int varid)
{
    dtf_var_t *ptr = varlist;

    while(ptr != NULL)
    {
        if(ptr->id == varid)
            break;

        ptr = ptr->next;
    }

    return ptr;
}

void add_var(dtf_var_t **vars, dtf_var_t *var)
{
    if(*vars == NULL)
        *vars = var;
    else{
        dtf_var_t *tmp = *vars;
        while(tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = var;
    }
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

static void delete_var(dtf_var_t *var)
{
   buffer_node_t *node = var->nodes;
    while(node != NULL){
        dtf_free(node->data, node->data_sz);
        var->nodes = var->nodes->next;
        dtf_free(node, sizeof(buffer_node_t));
        node = var->nodes;
    }
    dtf_free(var->shape, var->ndims*sizeof(MPI_Offset));
    if(var->distr_count != NULL)
        dtf_free(var->distr_count, var->ndims*sizeof(MPI_Offset));
    if(var->first_coord != NULL)
        dtf_free(var->first_coord, var->ndims*sizeof(MPI_Offset));
    dtf_free(var, sizeof(dtf_var_t));
}

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* fbuf)
{
    dtf_var_t    *var;
    if(fbuf == NULL)
        return;

    assert(buflist != NULL);

    if(fbuf->header != NULL)
        dtf_free(fbuf->header, fbuf->hdr_sz);

    var = fbuf->vars;
    while(var != NULL){
        fbuf->vars = var->next;
        delete_var(var);
        var = fbuf->vars;
    }

    assert(fbuf->ioreqs == NULL);
    assert(fbuf->rreq_cnt == 0);
    assert(fbuf->wreq_cnt == 0);
    if(fbuf->mst_info->iodb != NULL){
        clean_iodb(fbuf->mst_info->iodb);
        dtf_free(fbuf->mst_info->iodb, sizeof(ioreq_db_t));
    }

    dtf_free(fbuf->mst_info->masters, fbuf->mst_info->nmasters*sizeof(int));
    dtf_free(fbuf->mst_info, sizeof(master_info_t));
    dtf_free(fbuf->distr_ranks, fbuf->distr_nranks*sizeof(int));

//    if(*buflist == fbuf)
//        *buflist = fbuf->next;
//    else{
//        prev = *buflist;
//        while(prev->next != fbuf)
//            prev = prev->next;
//        prev->next = fbuf->next;
//    }

    dtf_free(fbuf, sizeof(file_buffer_t));
}

file_buffer_t* new_file_buffer()
{
    file_buffer_t *buf;

    buf = (file_buffer_t*)dtf_malloc(sizeof(struct file_buffer));
    assert( buf != NULL );
    buf->file_path[0]='\0';
    buf->alias_name[0]='\0';
    buf->next = NULL;
    buf->reader_id = -1;
    buf->version = 0;
    buf->writer_id=-1;
    buf->is_ready = 0;
    buf->vars = NULL;
    buf->var_cnt = 0;
    buf->header = NULL;
    buf->hdr_sz = 0;
    buf->iomode = DTF_UNDEFINED;
    buf->distr_rule = DISTR_RULE_P2P;
    buf->distr_pattern = DISTR_PATTERN_ALL;
    buf->distr_range = 0;
    buf->distr_nranks = 0;
    buf->distr_ranks = NULL;
    buf->distr_ndone = 0;
    //buf->ioreq_cnt = 0;
    buf->rreq_cnt = 0;
    buf->wreq_cnt = 0;
    buf->ioreqs = NULL;
    buf->ncid = -1;
    buf->explicit_match = 0;
    buf->done_matching_flag = 0;
    buf->rdr_closed_flag = 0;
    buf->fready_notify_flag = DTF_UNDEFINED;
    buf->comm = MPI_COMM_NULL;
    buf->root_writer = -1;
    buf->root_reader = -1;
    buf->nwriters = 0;
    buf->mst_info = dtf_malloc(sizeof(master_info_t));
    assert(buf->mst_info != NULL);
    buf->mst_info->masters = NULL;
    buf->mst_info->is_master_flag = 0;
    buf->mst_info->my_workgroup_sz = 0;
    buf->mst_info->nmasters = 0;
    buf->mst_info->iodb = NULL;
    buf->mst_info->nrranks_completed = 0;
    buf->mst_info->nwranks_completed = 0;
    buf->mst_info->nrranks_opened = 0;
    buf->mst_info->nwranks_opened = 0;
    buf->is_matching_flag = 0;
    return buf;
}

int has_unlim_dim(dtf_var_t *var)
{
    int ret = 0, i;
    for(i = 0; i < var->ndims; i++)
        if(var->shape[i] == DTF_UNLIMITED){
            ret = 1;
            break;
    }
    return ret;
}

dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i;
    dtf_var_t *var = (dtf_var_t*)dtf_malloc(sizeof(dtf_var_t));
    assert(var!=NULL);

    /*Initialize whatever we can initialize at this stage*/
    var->nodes = NULL;
    var->node_cnt=0;
    var->id = varid;
    if(ndims > 0){ //non scalar
        var->shape = (MPI_Offset*)dtf_malloc(ndims*sizeof(MPI_Offset));
        for(i=0; i<ndims;i++)
            var->shape[i] = shape[i];
    }
    else
        var->shape = NULL;
    var->first_coord = NULL;
    var->ndims = ndims;
    var->distr_count = NULL;
    var->next = NULL;
    var->dtype = dtype;
    return var;
}


int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start, const MPI_Offset *count )
{
    dtf_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "Did not find var id %d in file %s", varid, fbuf->file_path);
        return 1;
    }

    if(var->ndims > 0){
        int i;

        if(frt_indexing){
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[var->ndims-i-1]){
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[var->ndims-i-1]);
                    return 1;
                }
        } else {
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[i]){
                            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[i]);
                            return 1;
                }
        }
    }
//    else {
//        if( (start != NULL) || (count != NULL)){
//            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d is a scalar variable but trying to read an array", varid);
//            return 1;
//        }
//    }
    return 0;
}
