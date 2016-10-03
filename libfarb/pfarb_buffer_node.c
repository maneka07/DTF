#include "pfarb_buffer_node.h"
#include <assert.h>
#include <string.h>
#include "pfarb_common.h"

void print_nodes(buffer_node_t* nodes){
    buffer_node_t *tmp = nodes;

    while(tmp != NULL){
        FARB_DBG(VERBOSE_ALL_LEVEL, "node (%lld, %lld)", tmp->offset, tmp->data_sz);
        tmp = tmp->next;
    }
}

buffer_node_t* new_buffer_node(MPI_Offset offset, MPI_Offset data_sz, int malloc_flag)
{
    buffer_node_t *node = malloc(sizeof(buffer_node_t));
    assert(node != NULL);
    if(malloc_flag){
        node->data = malloc((size_t)data_sz);
        assert(node->data != NULL);
    } else
        node->data = NULL;

    node->data_sz = data_sz;
    node->offset = offset;
    node->next = NULL;
    node->prev = NULL;

    return node;
}

void insert_buffer_node(buffer_node_t** list, buffer_node_t* node)
{
    FARB_DBG(VERBOSE_ALL_LEVEL, "Insert node offt %d, datasz %d", (int)node->offset, (int)node->data_sz);
    if(*list == NULL){
        *list = node;
        return;
    }else{
        buffer_node_t *tmp = *list;
        while( tmp->next != NULL){
            if(tmp->offset > node->offset)
                break;
            tmp = tmp->next;
        }
        if(tmp->offset > node->offset){
            if(tmp->prev != NULL){
                tmp->prev->next = node;
                node->prev = tmp->prev;
            }
            tmp->prev = node;
            node->next = tmp;
            if(tmp == *list)
                *list = node;
        } else {
            node->prev = tmp;
            node->next = tmp->next;
            tmp->next = node;
            if(node->next != NULL)
                node->next->prev = node;
        }
    }
}

/*Check if process has written some data at this offset*/
int data_present_at_offt(buffer_node_t *list, MPI_Offset offset)
{
    buffer_node_t *tmp = list;
    int is_present = 0;
    while(tmp != NULL){
        if( (offset >= tmp->offset) && (offset < tmp->offset+tmp->data_sz)){
            is_present = 1;
            break;
        }
        tmp = tmp->next;
    }
    return is_present;
}
