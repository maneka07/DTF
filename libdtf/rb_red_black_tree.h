﻿/*
 © Copyright 2019 RIKEN Center for Computational Science, System Software
 *	Development Team All rights reserved.

 * The Data Transfer Framework (DTF) was designed to work with 
 * the PnetCDF library developed by Northwestern University and 
 * Argonne National Laboratory. The modified version of PnetCDF is 
 * provided with this distribution.
 * 
 * Access and use of this software shall impose the following obligations 
 * and understandings on the user. The user is granted the right, without 
 * any fee or cost, to use, copy, modify, alter, enhance and distribute 
 * this software, and any derivative works thereof, and its supporting 
 * documentation for any purpose whatsoever, provided that this entire 
 * notice appears in all copies of the software, derivative works and 
 * supporting documentation.  Further, RIKEN requests that the user credit 
 * RIKEN in any publications that result from the use of this software or 
 * in any product that includes this software.  The name RIKEN, however, 
 * may not be used in any advertising or publicity to endorse or promote 
 * any products or commercial entity unless specific written permission is 
 * obtained from RIKEN. The user also understands that RIKEN is not 
 * obligated to provide the user with any support, consulting, training or 
 * assistance of any kind with regard to the use, operation and 
 * performance of this software nor to provide the user with any updates, 
 * revisions, new versions or "bug fixes."
 * 
 * THIS SOFTWARE IS PROVIDED BY RIKEN "AS IS" AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN 
 * NO EVENT SHALL RIKEN BE LIABLE FOR ANY SPECIAL, INDIRECT OR 
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF 
 * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR 
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, 
 * USE OR PERFORMANCE OF THIS SOFTWARE.
 *
*/

#ifndef RB_RED_BLACK_TREE_H_INCLUDED
#define RB_RED_BLACK_TREE_H_INCLUDED

#ifdef DMALLOC
#include <dmalloc.h>
#endif
#include"rb_misc.h"
#include"rb_stack.h"
#include <mpi.h>

/*  CONVENTIONS:  All data structures for red-black trees have the prefix */
/*                "rb_" to prevent name conflicts. */
/*                                                                      */
/*                Function names: Each word in a function name begins with */
/*                a capital letter.  An example funcntion name is  */
/*                CreateRedTree(a,b,c). Furthermore, each function name */
/*                should begin with a capital letter to easily distinguish */
/*                them from variables. */
/*                                                                     */
/*                Variable names: Each word in a variable name begins with */
/*                a capital letter EXCEPT the first letter of the variable */
/*                name.  For example, int newLongInt.  Global variables have */
/*                names beginning with "g".  An example of a global */
/*                variable name is gNewtonsConstant. */

/* comment out the line below to remove all the debugging assertion */
/* checks from the compiled code.  */
#define DEBUG_ASSERT 1

typedef struct rb_red_blk_node {
  void* key;
  void* info;
  int red; /* if red=0 then the node is black */
  struct rb_red_blk_node* left;
  struct rb_red_blk_node* right;
  struct rb_red_blk_node* parent;
} rb_red_blk_node;


/* Compare(a,b) should return 1 if *a > *b, -1 if *a < *b, and 0 otherwise */
/* Destroy(a) takes a pointer to whatever key might be and frees it accordingly */
typedef struct rb_red_blk_tree {
  int (*Compare)(const void* a, const void* b);
  void (*DestroyKey)(void* a);
  void (*DestroyInfo)(void* a);
  void (*PrintKey)(const void* a);
  void (*PrintInfo)(void* a);
  /*  A sentinel is used for root and for nil.  These sentinels are */
  /*  created when RBTreeCreate is caled.  root->left should always */
  /*  point to the node which is the root of the tree.  nil points to a */
  /*  node which should always be black but has aribtrary children and */
  /*  parent and no key or info.  The point of using these sentinels is so */
  /*  that the root and nil nodes do not require special cases in the code */
  rb_red_blk_node* root;
  rb_red_blk_node* nil;
  /*DTF*/
  int cur_dim;
  unsigned long nnodes;
} rb_red_blk_tree;


/*DTF structs*/

typedef struct block{
    int rank;
    MPI_Offset *start;
    MPI_Offset *count;
}block_t;

typedef struct insert_info{
	int cur_dim;
	int ndims;
	block_t *blck;		
}insert_info;

typedef struct node_info{
	MPI_Offset	max_node_rcoord;
	MPI_Offset	max_subtr_rcoord;
	block_t    *blck;
	struct rb_red_blk_tree *next_dim_tree;   //tree sorting blocks in the dimension cur_dim+1
}node_info;

rb_red_blk_tree* RBTreeCreate(int  (*CompFunc)(const void*, const void*),
			     void (*DestFunc)(void*),
			     void (*InfoDestFunc)(void*),
			     void (*PrintFunc)(const void*),
			     void (*PrintInfo)(void*));
rb_red_blk_tree* RBTreeCreateBlocks(int  (*CompFunc)(const void*, const void*),
			     void (*DestFunc)(void*),
			     void (*InfoDestFunc)(void*),
			     void (*PrintFunc)(const void*),
			     void (*PrintInfo)(void*), int cur_dim);
rb_red_blk_node * RBTreeInsert(rb_red_blk_tree*, void* key, void* info);
rb_red_blk_node * RBTreeInsertBlock(rb_red_blk_tree* tree, void* ins_info);

void RBTreePrint(rb_red_blk_tree*);
void RBTreePrintBlocks(rb_red_blk_tree *tree);
void RBDelete(rb_red_blk_tree* , rb_red_blk_node* );
void RBTreeDestroy(rb_red_blk_tree*);
rb_red_blk_node* TreePredecessor(rb_red_blk_tree*,rb_red_blk_node*);
rb_red_blk_node* TreeSuccessor(rb_red_blk_tree*,rb_red_blk_node*);
rb_red_blk_node* RBExactQuery(rb_red_blk_tree*, void*);
rb_red_blk_node* RBExactQueryBlock(rb_red_blk_tree* tree, MPI_Offset* start, int ndims);
stk_stack * RBEnumerate(rb_red_blk_tree* tree,void* low, void* high);
void NullFunction(void*);

/***********************DTF***************************************/
void rb_print_info(void *ninfo);
void rb_print_key(const void *key);
int rb_key_cmp(const void *a, const void *b);
void rb_destroy_node_info(void* ninfo);

void rb_print_stats();
block_t *rb_find_block(rb_red_blk_tree* tree, MPI_Offset* start, MPI_Offset* count, int ndims);
#endif
