#include <assert.h>
#include <string.h>
#include "pfarb_req_match.h"
#include "pfarb_file_buffer.h"
#include "pfarb.h"
#include "pfarb_common.h"


file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid)
{
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
        if(file_path == NULL){
            if(ptr->ncid == ncid)
                break;
        } else {
            if(strstr(file_path, ptr->file_path) != NULL || ( strlen(ptr->alias_name) != 0 && strstr(file_path, ptr->alias_name) != NULL) )
                break;
        }
        ptr = ptr->next;
    }

    return ptr;
}

farb_var_t* find_var(farb_var_t* varlist, int varid)
{
    farb_var_t *ptr = varlist;

    while(ptr != NULL)
    {
        if(ptr->id == varid)
            break;

        ptr = ptr->next;
    }

    return ptr;
}

void add_var(farb_var_t **vars, farb_var_t *var)
{
    if(*vars == NULL)
        *vars = var;
    else{
        farb_var_t *tmp = *vars;
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

static void delete_var(farb_var_t *var)
{
   buffer_node_t *node = var->nodes;

    while(node != NULL){
        free(node->data);
        var->nodes = var->nodes->next;
        free(node);
        node = var->nodes;
    }

    if(var->distr_count != NULL)
        free(var->distr_count);
    if(var->first_coord != NULL)
        free(var->distr_count);
    free(var);
}

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* fbuf)
{

    file_buffer_t *prev;
    farb_var_t    *var, *tmp;
    if(fbuf == NULL)
        return;

    assert(buflist != NULL);

    if(fbuf->header != NULL)
        free(fbuf->header);

    var = fbuf->vars;
    while(var != NULL){
        tmp = var->next;
        delete_var(var);
        var = tmp;
    }

    if(fbuf->iodb != NULL){
        clean_iodb(fbuf->iodb);
        free(fbuf->iodb);
    }

    free(fbuf->distr_ranks);

    if(*buflist == fbuf)
        *buflist = fbuf->next;
    else{
        prev = *buflist;
        while(prev->next != fbuf)
            prev = prev->next;
        prev->next = fbuf->next;
    }

    free(fbuf);
}

file_buffer_t* new_file_buffer()
{
    file_buffer_t *buf;

    buf = (file_buffer_t*)malloc(sizeof(struct file_buffer));
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
    buf->iomode = FARB_UNDEFINED;
    buf->distr_rule = DISTR_RULE_P2P;
    buf->distr_pattern = DISTR_PATTERN_ALL;
    buf->distr_range = 0;
    buf->distr_nranks = 0;
    buf->distr_ranks = NULL;
    buf->distr_ndone = 0;
    buf->hdr_sent_flag = 0;
    buf->ioreq_cnt = 0;
    buf->ioreqs = NULL;
    buf->iodb = NULL;
    buf->ncid = -1;
    buf->explicit_match = 0;
    buf->done_matching_flag = 0;
    buf->fclosed_flag = 0;
    return buf;
}

int has_unlim_dim(farb_var_t *var)
{
    int ret = 0, i;
    for(i = 0; i < var->ndims; i++)
        if(var->shape[i] == FARB_UNLIMITED){
            ret = 1;
            break;
    }
    return ret;
}

farb_var_t* new_var(int varid, int ndims, MPI_Offset el_sz, MPI_Offset *shape)
{
    int i;
    farb_var_t *var = malloc(sizeof(farb_var_t));
    assert(var!=NULL);

    /*Initialize whatever we can initialize at this stage*/
    var->nodes = NULL;
    var->node_cnt=0;
    var->id = varid;
    if(ndims > 0){ //non scalar
        var->shape = (MPI_Offset*)malloc(ndims*sizeof(MPI_Offset));
        for(i=0; i<ndims;i++)
            var->shape[i] = shape[i];
    }
    else
        var->shape = NULL;
    var->first_coord = NULL;
    var->ndims = ndims;
    var->distr_count = NULL;
    var->next = NULL;
    var->el_sz = el_sz;
    return var;
}


int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start, const MPI_Offset *count )
{
    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "Did not find var id %d in file %s", varid, fbuf->file_path);
        return 1;
    }

    if(var->ndims > 0){
        int i;
        for(i = 0; i < var->ndims; i++)
            if(var->shape[i] == FARB_UNLIMITED) //no boundaries for unlimited dimension
                continue;
            else if(start[i] + count[i] > var->shape[i]){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[i]);
                return 1;
            }
    } else {
        if( (start != NULL) || (count != NULL)){
            FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: var %d is a scalar variable but trying to read an array", varid);
            return 1;
        }
    }
    return 0;
}
