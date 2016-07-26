#include "farb_file_buffer.h"
#include <assert.h>
file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path)
{
    struct file_buffer *ptr = buflist;
//    FARB_DBG(VERBOSE_DBG_LEVEL,   "Will search for %s", file_path);
    while(ptr != NULL)
    {
        if(strstr(file_path, ptr->file_path) != NULL || ( strlen(ptr->alias_name) != 0 && strstr(file_path, ptr->alias_name) != NULL) )
            break;

        ptr = ptr->next;
    }
//    if(ptr==NULL)
//        FARB_DBG(VERBOSE_DBG_LEVEL,   "File %s not found", file_path);
//    else
//        FARB_DBG(VERBOSE_DBG_LEVEL,   "File found %s", ptr->file_path);

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

void delete_file_buffer(file_buffer_t** buflist, file_buffer_t* buf)
{
    int i;
    file_buffer_t *prev;
    if(buf == NULL)
        return;

    assert(buflist != NULL);

    if(buf->reader_ids != NULL)
        free(buf->reader_ids);

    for(i = 0; i < buf->node_cnt; i++){
        free(buf->nodes_tbl[i].data);
    }
    if(buf->nodes_tbl != NULL)
        free(buf->nodes_tbl);
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
        buf->nodes = NULL;
        buf->data_sz = 0;
        buf->node_cnt = 0;
        buf->nodes_tbl = NULL;
        buf->mode = IO_MODE_UNDEFINED;
    } else {
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "Error allocating memory");
    }

    return buf;
}
