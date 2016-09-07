#include "pfarb_mem.h"
#include <assert.h>

static void print_nodes(buffer_node_t* nodes){
    buffer_node_t *tmp = nodes;

    while(tmp != NULL){
        FARB_DBG(VERBOSE_DBG_LEVEL, "node (%lld, %lld)", tmp->offset, tmp->data_sz);
        tmp = tmp->next;
    }
}

MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data)
{
    MPI_Offset offt, dsz, to_cpy, copied=0l;
    buffer_node_t *tmp;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Enter mem_write");
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Will write %llu bytes to var %d at offt %llu", data_sz, var->id, offset);

    if(var->nodes == NULL){
        buffer_node_t *node = new_buffer_node(offset, data_sz);
        var->node_cnt++;
        memcpy(node->data, data, (size_t)data_sz);
        insert_buffer_node(&(var->nodes), node);
        FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
        return data_sz;
    }

    print_nodes(var->nodes);

    tmp = var->nodes;

    while(tmp->next != NULL){
        if( (offset < tmp->offset) || (offset < tmp->offset+tmp->data_sz))
            break;
        tmp = tmp->next;
    }

    offt = offset;
    dsz = data_sz;

    while(copied != data_sz){

        FARB_DBG(VERBOSE_DBG_LEVEL, "tmp %p var %p", tmp, var->nodes);
        assert(tmp != NULL);

        if(offt < tmp->offset){
            FARB_DBG(VERBOSE_DBG_LEVEL, "1");
            to_cpy = tmp->offset - offt;
            if(to_cpy > dsz)
                to_cpy = dsz;
            buffer_node_t *node = new_buffer_node(offset, to_cpy);
            var->node_cnt++;
            memcpy(node->data, data+copied, (size_t)to_cpy);
            insert_buffer_node(&(var->nodes), node);

        } else if(offt >= tmp->offset + tmp->data_sz){ /*special case: the offset is > the last node*/
            FARB_DBG(VERBOSE_DBG_LEVEL, "2");
            buffer_node_t *node = new_buffer_node(offset, dsz);
            var->node_cnt++;
            to_cpy = dsz;
            memcpy(node->data, data+copied, (size_t)to_cpy);
            insert_buffer_node(&(var->nodes), node);

        } else { /*overlapping*/
            FARB_DBG(VERBOSE_DBG_LEVEL, "3");
            //do we need to extend the buffer?
            if(offt + dsz > tmp->offset + tmp->data_sz){
                MPI_Offset ext_sz = (offt + dsz) - (tmp->offset + tmp->data_sz);
                if( (tmp->next != NULL) && (offt+dsz > tmp->next->offset))
                    ext_sz = tmp->next->offset - (tmp->offset+tmp->data_sz); //adjust
                tmp->data = realloc(tmp->data, tmp->data_sz + ext_sz);
                tmp->data_sz += ext_sz;
                assert(tmp->data != NULL);
                to_cpy = tmp->offset + tmp->data_sz - offt + ext_sz;
            } else
                to_cpy = dsz;

            memcpy(tmp->data+(offt - tmp->offset), data+copied, (size_t)to_cpy);
        }

        copied += to_cpy;
        FARB_DBG(VERBOSE_DBG_LEVEL, "Copied %d, left to copy %d", (int)copied, (int)(data_sz - copied));
        offt += to_cpy;
        dsz -= to_cpy;
        tmp = tmp->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
    return copied;
}

MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data)
{

    MPI_Offset copied, node_offt, to_cpy;


    FARB_DBG(VERBOSE_DBG_LEVEL, "mem_read at offt %ld of size %lu", (long int)offset, (long unsigned)data_sz);

    //Find the first node that holds the data
    buffer_node_t *tmp = var->nodes;

    print_nodes(var->nodes);

    while(tmp->next != NULL){
        if( (offset >= tmp->offset) && (offset < tmp->offset+tmp->data_sz))
            break;
        tmp = tmp->next;
    }

    if(tmp == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB warning: there is no data at offset %llu", offset);
        return 0;
    }
    copied = 0;

    while(copied != data_sz){
        node_offt = offset - tmp->offset;

        if(tmp->data_sz - node_offt >= data_sz)
            to_cpy = data_sz;
        else
            to_cpy = tmp->data_sz - node_offt;

        memcpy(data+copied, tmp->data+node_offt, (size_t)to_cpy);
        copied += to_cpy;

        tmp = tmp->next;
        if(tmp == NULL)
            break;

        if( (offset+copied < tmp->offset) || (offset+copied >= tmp->offset+tmp->data_sz) )
            break;
    }

    return copied;
}
