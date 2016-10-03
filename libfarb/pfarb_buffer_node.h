#ifndef BUFFER_NODE_H_INCLUDED
#define BUFFER_NODE_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define DEFAULT_BUFFER_NODE_SIZE (1024*1024)


typedef struct buffer_node{
    MPI_Offset offset;      /*offset in bytes from the beginning of the variable*/
    MPI_Offset data_sz;
    void* data;
    struct buffer_node *next;
    struct buffer_node *prev;
}buffer_node_t;

buffer_node_t* new_buffer_node(MPI_Offset offset, MPI_Offset data_sz, int malloc_flag);
buffer_node_t* find_buffer_node(buffer_node_t* nodelist, MPI_Offset offset);
void insert_buffer_node(buffer_node_t** list, buffer_node_t* node);
int data_present_at_offt(buffer_node_t *list, MPI_Offset offset);
void print_nodes(buffer_node_t* nodes);
#endif /*BUFFER_NODE_H_INCLUDED*/
