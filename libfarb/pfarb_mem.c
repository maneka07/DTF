#include <assert.h>
#include <string.h>
#include "pfarb_mem.h"
#include "pfarb.h"
#include "pfarb_common.h"
#include "pfarb_util.h"


//static MPI_Offset mem_recursive_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data, int dim, MPI_Offset coord[])
//{
//    int i;
//    MPI_Offset written = 0;
//    if(dim == var->ndims - 1){
//        MPI_Offset *dst_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
//        assert(dst_coord != NULL);
//        MPI_Offset data_sz = count[var->ndims-1]*var->el_sz;
//        //printf("%d %d %d\n", coord[0], coord[1], coord[2]);
//        MPI_Offset src_start_idx = to_1d_index(var->ndims, count, coord); //el offset within *data
//        MPI_Offset src_offt = src_start_idx * var->el_sz;
//        for(i = 0; i < var->ndims; i++)
//          dst_coord[i] = coord[i] + start[i];
//        FARB_DBG(VERBOSE_ALL_LEVEL, "will write %llu bytes from el %llu (%llu)", data_sz, src_start_idx, src_offt);
//        written = mem_write(var, dst_coord, data_sz, data+src_offt);
//        free(dst_coord);
//        return written;
//    }
//
//    for(i = 0; i < count[dim]; i++){
//        coord[dim] = i;
//        written += mem_recursive_write(var, start, count, data, dim+1, coord);
//    }
//    return written;
//}

MPI_Offset mem_noncontig_write(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data)
{
    MPI_Offset written = 0;
    MPI_Offset src_offset = 0;
    contig_mem_chunk_t *chunks, *tmp;
    int el_sz;

    MPI_Type_size(var->dtype, &el_sz);

    int nelems = 0;
    get_contig_mem_list(var, var->dtype, start, count, &nelems, &chunks);

    //for(i = 0; i < nelems; i = i+2){
    while(chunks != NULL){
        written += mem_contiguous_write(var, chunks->offset, chunks->data_sz, (unsigned char*)data+src_offset);
        src_offset += chunks->data_sz;
        tmp = chunks;
        chunks = chunks->next;
        free(tmp);
    }
    if((var->ndims > 0) && (nelems > 0)){
        if(var->first_coord == NULL){
            var->first_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
            assert(var->first_coord != NULL);
            memcpy((void*)var->first_coord, (void*)start, var->ndims*sizeof(MPI_Offset));
        } else {
            MPI_Offset idx1 = to_1d_index(var->ndims, var->shape, var->first_coord);
            MPI_Offset idx2 = to_1d_index(var->ndims, var->shape, start);
            if(idx2 < idx1){
                memcpy((void*)var->first_coord, (void*)start, var->ndims*sizeof(MPI_Offset));
                assert(var->nodes->offset == idx2*el_sz);
            }
        }
    }

    return written;
}

MPI_Offset mem_contiguous_write(farb_var_t *var, MPI_Offset offset, MPI_Offset data_sz, void *data)
{
    MPI_Offset left, to_cpy, copied=0;
    buffer_node_t *node;
    FARB_DBG(VERBOSE_ALL_LEVEL,   "Will write %llu bytes to offt %llu", data_sz, offset);

    if(var->nodes == NULL){
        buffer_node_t *new_node = new_buffer_node(offset, data_sz, 1);
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
            buffer_node_t *new_node = new_buffer_node(offset, to_cpy, 1);
            var->node_cnt++;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt 0",
                     (int)to_cpy, (int)new_node->offset, (int)new_node->data_sz);
            memcpy(new_node->data, (unsigned char*)data+copied, (size_t)to_cpy);
            insert_buffer_node(&(var->nodes), new_node);

        } else if(offset > node->offset + node->data_sz){ /*the offset is > the last node*/
            buffer_node_t *new_node = new_buffer_node(offset, left, 1);
            var->node_cnt++;
            to_cpy = left;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Copy %d to node with offt %d (sz %d) at node offt 0",
                     (int)to_cpy, (int)new_node->offset, (int)new_node->data_sz);
            memcpy(new_node->data, (unsigned char*)data+copied, (size_t)to_cpy);
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

            memcpy( (unsigned char*)node->data+node->data_sz,  (unsigned char*)data+copied, (size_t)to_cpy);
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
            memcpy( (unsigned char*)node->data+(offset - node->offset),  (unsigned char*)data+copied, (size_t)to_cpy);
        }

        copied += to_cpy;
        offset += to_cpy;
        left -= to_cpy;
        node = node->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
    return copied;
}

static void traverse_dims(farb_var_t *var,
                          MPI_Datatype dtype,
                          const MPI_Offset count[],
                          const MPI_Offset start[],
                          int dim,
                          MPI_Offset coord[],
                          int *nelems,
                          contig_mem_chunk_t **list)
{
    int i;
    if(dim == var->ndims - 1){
        int def_el_sz, usr_el_sz;
        contig_mem_chunk_t *new_chunk = (contig_mem_chunk_t*)malloc(sizeof(contig_mem_chunk_t));
        assert(new_chunk!=NULL);

        if(*list == NULL){
            *list = new_chunk;
        }
        else{
            contig_mem_chunk_t *tmp = *list;
            while(tmp->next != NULL)
                tmp = tmp->next;
            tmp->next = new_chunk;
        }
        MPI_Offset *tmp = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
        assert(tmp != NULL);
        for(i = 0; i < var->ndims; i++){
            tmp[i] = coord[i] - start[i];
            //FARB_DBG(VERBOSE_ALL_LEVEL, "coord %d", (int)coord[i]);
        }
        MPI_Type_size(var->dtype, &def_el_sz);
        MPI_Type_size(dtype, &usr_el_sz);
        new_chunk->offset = to_1d_index(var->ndims, var->shape, coord)*def_el_sz; //offset
        new_chunk->usrbuf_offset = to_1d_index(var->ndims, count, tmp)*usr_el_sz; //offset inside the user buffer
        new_chunk->data_sz = count[var->ndims-1]*def_el_sz;    //data_sz
        new_chunk->next = NULL;
        (*nelems)++;
        free(tmp);
        FARB_DBG(VERBOSE_ALL_LEVEL,"Added a tuple (off %d, uoff %d, sz %d)", (int)new_chunk->offset, (int)new_chunk->usrbuf_offset, (int)new_chunk->data_sz);
        return;
    }

    for(i = 0; i < count[dim]; i++){
        coord[dim] = start[dim] + i;
        traverse_dims(var, dtype, count, start, dim+1, coord, nelems, list);
    }
}

/*Returns an array consisting of tuples (offset, data_sz) of contiguous memory chunks that hold
  the block of data of a multi-dimensional var when it's flattened to a 1d array.
  The data block starts at coordinate start[] and is of size count[].*/
void get_contig_mem_list(farb_var_t *var,
                         MPI_Datatype dtype,
                         const MPI_Offset start[],
                         const MPI_Offset count[],
                         int *nelems,
                         contig_mem_chunk_t **list)
{
    MPI_Offset *cur_coord = (MPI_Offset*)malloc(var->ndims*sizeof(MPI_Offset));
    memcpy((void*)cur_coord, (void*)start, var->ndims*sizeof(MPI_Offset));
    *list = NULL;
    *nelems = 0;
    traverse_dims(var, dtype, count, start, 0, cur_coord, nelems, list);
    free(cur_coord);
}

MPI_Offset mem_noncontig_read(farb_var_t *var, const MPI_Offset start[], const MPI_Offset count[], void *data)
{
    MPI_Offset readsz = 0, dst_offset = 0;
    contig_mem_chunk_t *chunks, *tmp;
    int nelems = 0;
    get_contig_mem_list(var, var->dtype, start, count, &nelems, &chunks);

    while(chunks != NULL){
        readsz += mem_contiguous_read(var, chunks->offset, chunks->data_sz,  (unsigned char*)data+dst_offset);
        dst_offset += chunks->data_sz;
        tmp = chunks;
        chunks = chunks->next;
        free(tmp);
    }

    return readsz;
}

MPI_Offset mem_contiguous_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data)
{

    MPI_Offset copied, node_offt, to_cpy;
   // MPI_Offset start_idx = to_1d_index(var->ndims, var->shape, first_el_coord); //el offset within var->nodes
    //MPI_Offset offset = start_idx * var->el_sz;
    FARB_DBG(VERBOSE_ALL_LEVEL, "will read %llu bytes from offt %llu", data_sz, offset);

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
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: there is no data at offset %llu", offset);
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
        memcpy( (unsigned char*)data+copied,  (unsigned char*)tmp->data+node_offt, (size_t)to_cpy);
        copied += to_cpy;

        tmp = tmp->next;
        if(tmp == NULL)
            break;

        if( (offset+copied < tmp->offset) || (offset+copied >= tmp->offset+tmp->data_sz) )
            break;
    }

    return copied;
}
