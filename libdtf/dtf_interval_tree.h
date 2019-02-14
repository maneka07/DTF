#ifndef _DTF_INTERVAL_TREE_H_
#define _DTF_INTERVAL_TREE_H_

#include <mpi.h>

struct interval_tree;

typedef struct interval_node {
  int id;
  MPI_Offset lkey;  /*start value of interval*/
  MPI_Offset rkey;  /*end value of interval*/
  struct interval_tree *next_dim_tree;
  struct block_ref* bl_refs;
  int height;  /*tree height*/
  struct interval_node* left;
  struct interval_node* right;
  struct interval_node* parent;
} interval_node_t;


typedef struct interval_tree{
  struct interval_node* root;
  int cur_dim;       /*for k-dimensional interval tree, what dimension is this tree for*/
  unsigned long nnodes;
} interval_tree_t;

/*list of written blocks and their ranks*/
typedef struct block{
	int id;
    int rank;
    MPI_Offset *start;
    MPI_Offset *count;
    struct block* next;
}block_t;

typedef struct insert_info{
  int ndims;
  struct block* bl;
}insert_info_t;

typedef struct block_ref{
	struct block* bl;
	struct block_ref* next;
}block_ref_t;

interval_tree_t* 	IntervalTreeCreate();
void 				IntervalTreeDestroy(interval_tree_t* tree);
void 				IntervalTreeInsert(interval_tree_t* tree, insert_info_t* info);
void 				IntervalTreePrint(interval_tree_t* tree);
block_t* 			IntervalTreeFindOverlap(interval_tree_t* tree, MPI_Offset* start, MPI_Offset *count, int ndims);
void 				IntervalTreePrintMem();
#endif
