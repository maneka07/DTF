#include "pfarb_mem.h"
#include <assert.h>

static void print_nodes(buffer_node_t* nodes){
    buffer_node_t *tmp = nodes;

    while(tmp != NULL){
        FARB_DBG(VERBOSE_ALL_LEVEL, "node (%lld, %lld)", tmp->offset, tmp->data_sz);
        tmp = tmp->next;
    }
}

static MPI_Offset mem_recursive_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data, int dim, MPI_Offset coord[])
{
    int i;
    MPI_Offset written = 0;
    if(dim == var->ndims - 1){
        MPI_Offset *dst_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
        assert(dst_coord != NULL);
        MPI_Offset data_sz = count[var->ndims-1]*var->el_sz;
        //printf("%d %d %d\n", coord[0], coord[1], coord[2]);
        MPI_Offset src_start_idx = to_1d_index(var->ndims, count, coord); //el offset within *data
        MPI_Offset src_offt = src_start_idx * var->el_sz;
        for(i = 0; i < var->ndims; i++)
          dst_coord[i] = coord[i] + start[i];
        FARB_DBG(VERBOSE_ALL_LEVEL, "will write %llu bytes from el %llu (%llu)", data_sz, src_start_idx, src_offt);
        written = mem_write(var, dst_coord, data_sz, data+src_offt);
        free(dst_coord);
        return written;
    }

    for(i = 0; i < count[dim]; i++){
        coord[dim] = i;
        written += mem_recursive_write(var, start, count, data, dim+1, coord);
    }
    return written;
}

MPI_Offset mem_noncontig_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data)
{
    MPI_Offset written;
    MPI_Offset *start_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
    memset((void*)start_coord, 0, var->ndims*sizeof(MPI_Offset));
    written = mem_recursive_write(var, start, count, data, 0, start_coord);
    free(start_coord);
    return written;
}

MPI_Offset mem_write(farb_var_t *var, MPI_Offset first_el_coord[],  MPI_Offset data_sz, void *data)
{
    MPI_Offset left, to_cpy, copied=0;
    buffer_node_t *node;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Enter mem_write");
    MPI_Offset start_idx = to_1d_index(var->ndims, var->shape, first_el_coord); //el offset within var->nodes
    MPI_Offset offset = start_idx*var->el_sz;
    FARB_DBG(VERBOSE_ALL_LEVEL,   "Will write %llu bytes to el %llu (%llu)", data_sz, start_idx, offset);

    if(var->nodes == NULL){
        buffer_node_t *new_node = new_buffer_node(offset, data_sz, var->ndims, first_el_coord, 1);
        var->node_cnt++;
        memcpy(new_node->data, data, (size_t)data_sz);
        insert_buffer_node(&(var->nodes), new_node);
        FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
        return data_sz;
    }

    print_nodes(var->nodes);

    node = var->nodes;

    while(node->next != NULL){
        if( (offset < node->offset) || (offset <= node->offset + node->data_sz))
            break;
        node = node->next;
    }

    left = data_sz;

    FARB_DBG(VERBOSE_ALL_LEVEL, "found node with offt %d", (int)node->offset);

    while(copied != data_sz){

        assert(node != NULL);

        if(offset < node->offset){
            to_cpy = node->offset - offset;
            if(to_cpy > left)
                to_cpy = left;
            buffer_node_t *new_node = new_buffer_node(offset, to_cpy, var->ndims, first_el_coord, 1);
            var->node_cnt++;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt 0",
                     (int)to_cpy, (int)new_node->offset, (int)new_node->data_sz);
            memcpy(new_node->data, data+copied, (size_t)to_cpy);
            insert_buffer_node(&(var->nodes), new_node);

        } else if(offset > node->offset + node->data_sz){ /*the offset is > the last node*/
            buffer_node_t *new_node = new_buffer_node(offset, left, var->ndims, first_el_coord, 1);
            var->node_cnt++;
            to_cpy = left;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt 0",
                     (int)to_cpy, (int)new_node->offset, (int)new_node->data_sz);
            memcpy(new_node->data, data+copied, (size_t)to_cpy);
            insert_buffer_node(&(var->nodes), new_node);
        } else if(offset == node->offset + node->data_sz){
            //need to extend the buffer
            void *tmp;
            MPI_Offset ext_sz = left;

            if( (node->next != NULL) && (offset + left > node->next->offset)){
               ext_sz = node->next->offset - offset;
            }
            tmp = realloc(node->data, node->data_sz + ext_sz);
            assert(tmp != NULL);
            node->data = tmp;
            to_cpy = ext_sz;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt %d",
                    (int)to_cpy, (int)node->offset, (int)node->data_sz+(int)ext_sz, (int)(offset - node->offset));

            memcpy(node->data+node->data_sz, data+copied, (size_t)to_cpy);
            node->data_sz += ext_sz;
        } else { /*overlapping*/
            //do we need to extend the buffer?
            if(offset + left > node->offset + node->data_sz){
                void *tmp;
                MPI_Offset ext_sz = (offset + left) - (node->offset + node->data_sz);
                if( (node->next != NULL) && (offset+left > node->next->offset))
                    ext_sz = node->next->offset - (node->offset + node->data_sz); //adjust
                tmp = realloc(node->data, node->data_sz + ext_sz);
                assert(tmp != NULL);
                node->data = tmp;
                node->data_sz += ext_sz;
                to_cpy = node->offset + node->data_sz - offset + ext_sz;
            } else
                to_cpy = left;

            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt %d",
                    (int)to_cpy, (int)node->offset, (int)node->data_sz, (int)(offset - node->offset));
            memcpy(node->data+(offset - node->offset), data+copied, (size_t)to_cpy);
        }

        copied += to_cpy;
        offset += to_cpy;
        left -= to_cpy;
        node = node->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
    return copied;
}



static MPI_Offset mem_recursive_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data, int dim, MPI_Offset coord[])
{
    int i;
    MPI_Offset read = 0;
    if(dim == var->ndims - 1){
        MPI_Offset *src_coord = malloc(var->ndims*sizeof(MPI_Offset));
        assert(src_coord != NULL);
        MPI_Offset data_sz = count[var->ndims-1]*var->el_sz;
        for(i = 0; i < var->ndims; i++){
          src_coord[i] = coord[i] + start[i];
          FARB_DBG(VERBOSE_ALL_LEVEL, "coord %d, start %d", (int)coord[i], (int)start[i]);
        }
        MPI_Offset dst_start_idx = to_1d_index(var->ndims, count, coord); //el offset within *data
        MPI_Offset dst_offt = dst_start_idx*var->el_sz;

        FARB_DBG(VERBOSE_ALL_LEVEL, "will read %llu bytes to dst el %d(%d)", data_sz, (int)dst_start_idx, (int)dst_offt);
        read = mem_read(var, src_coord, data_sz, data+dst_offt);
        free(src_coord);
        return read;
    }

    for(i = 0; i < count[dim]; i++){
        coord[dim] = i;
        read += mem_recursive_read(var, start, count, data, dim+1, coord);
    }
    return read;
}
MPI_Offset mem_noncontig_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data)
{
    MPI_Offset readsz;
    MPI_Offset *cur_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
    memset((void*)cur_coord, 0, var->ndims*sizeof(MPI_Offset));
    readsz = mem_recursive_read(var, start, count, data, 0, cur_coord);
    free(cur_coord);
    return readsz;
}

MPI_Offset mem_read(farb_var_t *var, MPI_Offset first_el_coord[],  MPI_Offset data_sz, void *data)
{

    MPI_Offset copied, node_offt, to_cpy;
    MPI_Offset start_idx = to_1d_index(var->ndims, var->shape, first_el_coord); //el offset within var->nodes
    MPI_Offset offset = start_idx * var->el_sz;
    FARB_DBG(VERBOSE_ALL_LEVEL, "will read %llu bytes from src el %llu (%llu)", data_sz, start_idx, offset);

    //Find the first node that holds the data
    buffer_node_t *tmp = var->nodes;
    assert(tmp != NULL);
    print_nodes(var->nodes);

    while(tmp != NULL){
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
        node_offt = offset + copied - tmp->offset;

        if(tmp->data_sz - node_offt >= data_sz)
            to_cpy = data_sz;
        else
            to_cpy = tmp->data_sz - node_offt;

        FARB_DBG(VERBOSE_ALL_LEVEL, "From node with offt %d Will copy %d bytes from node_offt %d", (int)tmp->offset, (int)to_cpy, (int)node_offt);
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
