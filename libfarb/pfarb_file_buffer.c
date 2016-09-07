#include "pfarb_file_buffer.h"
#include "pfarb.h"
#include <assert.h>
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path)
{
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
        if(strstr(file_path, ptr->file_path) != NULL || ( strlen(ptr->alias_name) != 0 && strstr(file_path, ptr->alias_name) != NULL) )
            break;

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
    free(var);
}

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf)
{

    file_buffer_t *prev;
    farb_var_t    *var, *tmp;
    if(buf == NULL)
        return;

    assert(buflist != NULL);

    if(buf->header != NULL)
        free(buf->header);

    var = buf->vars;
    while(var != NULL){
        tmp = var->next;
        delete_var(var);
        var = tmp;
    }

    if(*buflist == buf)
        *buflist = buf->next;
    else{
        prev = *buflist;
        while(prev->next != buf)
            prev = prev->next;
        prev->next = buf->next;
    }

    free(buf);
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
    buf->mode = FARB_UNDEFINED;
    buf->distr_pattern = FARB_UNDEFINED;
    buf->distr_range = 0;
    buf->distr_nranks = 0;
    buf->distr_ndone = 0;
   // buf->distr_ranks_expr = NULL;
    buf->distr_rule = DISTR_RULE_P2P;
    return buf;
}

farb_var_t* new_var(int varid, int ndims, MPI_Offset *shape)
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

    var->ndims = ndims;
    var->distr_count = NULL;
    var->next = NULL;
    //TODO !!!! when dims are registered define xsz size
//    var->xsz = 0;
    return var;
}


