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
    long long unsigned i;

    for(i = 0; i< var->bufnode_cnt; i++){
        if(var->bufnodes[i].data != NULL)
            free(var->bufnodes[i].data);
    }

    free(var);
}

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf)
{

    file_buffer_t *prev;
    farb_var_t    *var, *tmp;
    if(buf == NULL)
        return;

    assert(buflist != NULL);

    if(buf->reader_ids != NULL)
        free(buf->reader_ids);

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

    buf = (file_buffer_t*) malloc(sizeof(struct file_buffer));

    if(buf != NULL){
        buf->file_path[0]='\0';
        buf->alias_name[0]='\0';
        buf->next = NULL;
        buf->reader_ids = NULL;
        buf->version = 0;
        buf->writer_id=-1;
        buf->nreaders = 0;
        buf->is_ready = 0;
        buf->transfered = 0;
        buf->vars = NULL;
        buf->header = NULL;
        buf->hdr_sz = 0;
        buf->mode = FARB_IO_MODE_UNDEFINED;
    } else {
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "Error allocating memory");
    }

    return buf;
}

farb_var_t* new_var(int varid, int ndims, MPI_Offset *shape)
{
    int i;
    farb_var_t *var = malloc(sizeof(farb_var_t));
    assert(var!=NULL);

    /*Initialize whatever we can initialize at this stage*/
    var->bufnodes = NULL;
    var->bufnode_cnt = 0;
    var->id = varid;
    var->begin = -1;
    if(ndims > 0){ //non scalar
        var->shape = (MPI_Offset*)malloc(ndims*sizeof(MPI_Offset));
        for(i=0; i<ndims;i++)
            var->shape[i] = shape[i];
    }
    else
        var->shape = NULL;

    var->ndims = ndims;
    var->next = NULL;
    var->type = MPI_DATATYPE_NULL;
    var->xsz = 0;
    return var;
}


