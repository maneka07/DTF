#include "pfarb_buffer_node.h"
#include <assert.h>

buffer_node_t* new_buffer_node(MPI_Offset offset, MPI_Offset data_sz)
{
    buffer_node_t *node = malloc(sizeof(buffer_node_t));
    assert(node != NULL);
    node->data = malloc((size_t)data_sz);
    assert(node->data != NULL);
    node->data_sz = data_sz;
    node->next = NULL;
    node->offset = offset;
    node->prev = NULL;

    return node;
}

void insert_buffer_node(buffer_node_t** list, buffer_node_t* node)
{

    if(*list == NULL){
        *list = node;
        return;
    }else{
        buffer_node_t *tmp = *list;
        while( tmp->next != NULL){
            if(tmp->next->offset > node->offset)
                break;
            tmp = tmp->next;
        }

        node->prev = tmp;
        node->next = tmp->next;
        tmp->next = node;
        if(node->next != NULL)
            node->next->prev = node;
    }
}
