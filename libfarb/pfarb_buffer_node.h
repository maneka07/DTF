#ifndef BUFFER_NODE_H_INCLUDED
#define BUFFER_NODE_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

#define DEFAULT_BUFFER_NODE_SIZE (1024*1024)


typedef struct buffer_node{

    char* data;
    size_t data_sz;
    //TODO remove all the rest
	int offset_start;
	char *node_ptr;
  	int dirty_flag;
	struct buffer_node *next;
}buffer_node_t;

#endif /*BUFFER_NODE_H_INCLUDED
*/
