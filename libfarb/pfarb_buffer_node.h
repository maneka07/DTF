#ifndef BUFFER_NODE_H_INCLUDED
#define BUFFER_NODE_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

#define DEFAULT_BUFFER_NODE_SIZE (1024*1024)


typedef struct buffer_node{

    void* data;
    MPI_Offset data_sz;
}buffer_node_t;

#endif /*BUFFER_NODE_H_INCLUDED
*/
