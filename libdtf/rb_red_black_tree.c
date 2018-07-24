#include "rb_red_black_tree.h"
#include <assert.h>
#include "dtf_util.h"
rb_red_blk_node* FindBlockPoint(rb_red_blk_tree* tree, rb_red_blk_node *subtr, MPI_Offset* start, int ndims);
rb_red_blk_node* FindBlockOverlap(rb_red_blk_tree* tree, rb_red_blk_node *subtr, MPI_Offset* start, MPI_Offset* count, int ndims);

/***********************************************************************/
/*  FUNCTION:  RBTreeCreate */
/**/
/*  INPUTS:  All the inputs are names of functions.  CompFunc takes to */
/*  void pointers to keys and returns 1 if the first arguement is */
/*  "greater than" the second.   DestFunc takes a pointer to a key and */
/*  destroys it in the appropriate manner when the node containing that */
/*  key is deleted.  InfoDestFunc is similiar to DestFunc except it */
/*  recieves a pointer to the info of a node and destroys it. */
/*  PrintFunc recieves a pointer to the key of a node and prints it. */
/*  PrintInfo recieves a pointer to the info of a node and prints it. */
/*  If RBTreePrint is never called the print functions don't have to be */
/*  defined and NullFunction can be used.  */
/**/
/*  OUTPUT:  This function returns a pointer to the newly created */
/*  red-black tree. */
/**/
/*  Modifies Input: none */
/***********************************************************************/

unsigned long long nnodes=0;
unsigned long long ntrees=0;

rb_red_blk_tree* RBTreeCreate( int (*CompFunc) (const void*,const void*),
			      void (*DestFunc) (void*),
			      void (*InfoDestFunc) (void*),
			      void (*PrintFunc) (const void*),
			      void (*PrintInfo)(void*)) {
  rb_red_blk_tree* newTree;
  rb_red_blk_node* temp;

  newTree=(rb_red_blk_tree*) SafeMalloc(sizeof(rb_red_blk_tree));
  newTree->Compare=  CompFunc;
  newTree->DestroyKey= DestFunc;
  newTree->PrintKey= PrintFunc;
  newTree->PrintInfo= PrintInfo;
  newTree->DestroyInfo= InfoDestFunc;

  /*  see the comment in the rb_red_blk_tree structure in red_black_tree.h */
  /*  for information on nil and root */
  temp=newTree->nil= (rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  temp->parent=temp->left=temp->right=temp;
  temp->red=0;
  temp->key=NULL;
  
  temp=newTree->root=(rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  temp->parent=temp->left=temp->right=newTree->nil;
  temp->key=NULL;
  temp->red=0;
  return(newTree);
}

rb_red_blk_tree* RBTreeCreateBlocks( int (*CompFunc) (const void*,const void*),
			      void (*DestFunc) (void*),
			      void (*InfoDestFunc) (void*),
			      void (*PrintFunc) (const void*),
			      void (*PrintInfo)(void*), int cur_dim) {
  rb_red_blk_tree* newTree;
  rb_red_blk_node* temp;

  newTree=(rb_red_blk_tree*) SafeMalloc(sizeof(rb_red_blk_tree));
  ntrees++;
  newTree->Compare=  CompFunc;
  newTree->DestroyKey= DestFunc;
  newTree->PrintKey= PrintFunc;
  newTree->PrintInfo= PrintInfo;
  newTree->DestroyInfo= InfoDestFunc;
  
  newTree->cur_dim = cur_dim;
  newTree->nnodes = 0;
  /*  see the comment in the rb_red_blk_tree structure in red_black_tree.h */
  /*  for information on nil and root */
  temp=newTree->nil= (rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  temp->parent=temp->left=temp->right=temp;
  temp->red=0;
  temp->key=NULL;
  nnodes++;
  newTree->nnodes++;
  
  temp=newTree->root=(rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  temp->parent = temp->left=temp->right=newTree->nil;
  temp->key=NULL;
  temp->red=0;
  newTree->nnodes++;
  nnodes++;
  return(newTree);
}

/***********************************************************************/
/*  FUNCTION:  LeftRotate */
/**/
/*  INPUTS:  This takes a tree so that it can access the appropriate */
/*           root and nil pointers, and the node to rotate on. */
/**/
/*  OUTPUT:  None */
/**/
/*  Modifies Input: tree, x */
/**/
/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
/*            makes the parent of x be to the left of x, x the parent of */
/*            its parent before the rotation and fixes other pointers */
/*            accordingly. */
/***********************************************************************/

void LeftRotate(rb_red_blk_tree* tree, rb_red_blk_node* x) {
  MPI_Offset max_rcoord;
  rb_red_blk_node* y;
  rb_red_blk_node* nil=tree->nil;
	DTF_DBG(VERBOSE_RB_TREE_LEVEL,"rotate left at %lld\n", *(MPI_Offset*)x->key);
  /*  I originally wrote this function to use the sentinel for */
  /*  nil to avoid checking for nil.  However this introduces a */
  /*  very subtle bug because sometimes this function modifies */
  /*  the parent pointer of nil.  This can be a problem if a */
  /*  function which calls LeftRotate also uses the nil sentinel */
  /*  and expects the nil sentinel's parent pointer to be unchanged */
  /*  after calling this function.  For example, when RBDeleteFixUP */
  /*  calls LeftRotate it expects the parent pointer of nil to be */
  /*  unchanged. */

  y=x->right;
  x->right=y->left;

  if (y->left != nil) y->left->parent=x; /* used to use sentinel here */
  /* and do an unconditional assignment instead of testing for nil */

  y->parent=x->parent;

  /* instead of checking if x->parent is the root as in the book, we */
  /* count on the root sentinel to implicitly take care of this case */
  if( x == x->parent->left) {
    x->parent->left=y;
  } else {
    x->parent->right=y;
  }
  y->left=x;
  x->parent=y;

  /*Update rcoords*/
  max_rcoord = ((node_info*)(x->info))->max_node_rcoord;
  if( (x->left != nil) && (((node_info*)x->left->info)->max_subtr_rcoord > max_rcoord))
	max_rcoord = ((node_info*)x->left->info)->max_subtr_rcoord;
  if( (x->right != nil) && (((node_info*)x->right->info)->max_subtr_rcoord > max_rcoord))
	max_rcoord = ((node_info*)x->right->info)->max_subtr_rcoord;
  ((node_info*)x->info)->max_subtr_rcoord = max_rcoord;

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not red in LeftRotate");
#endif
}


/***********************************************************************/
/*  FUNCTION:  RighttRotate */
/**/
/*  INPUTS:  This takes a tree so that it can access the appropriate */
/*           root and nil pointers, and the node to rotate on. */
/**/
/*  OUTPUT:  None */
/**/
/*  Modifies Input?: tree, y */
/**/
/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
/*            makes the parent of x be to the left of x, x the parent of */
/*            its parent before the rotation and fixes other pointers */
/*            accordingly. */
/***********************************************************************/

void RightRotate(rb_red_blk_tree* tree, rb_red_blk_node* y) {
  rb_red_blk_node* x;
  rb_red_blk_node* nil=tree->nil;
  MPI_Offset max_rcoord;
  
  DTF_DBG(VERBOSE_RB_TREE_LEVEL,"rotate right at %lld\n", *(MPI_Offset*)y->key);
  /*  I originally wrote this function to use the sentinel for */
  /*  nil to avoid checking for nil.  However this introduces a */
  /*  very subtle bug because sometimes this function modifies */
  /*  the parent pointer of nil.  This can be a problem if a */
  /*  function which calls LeftRotate also uses the nil sentinel */
  /*  and expects the nil sentinel's parent pointer to be unchanged */
  /*  after calling this function.  For example, when RBDeleteFixUP */
  /*  calls LeftRotate it expects the parent pointer of nil to be */
  /*  unchanged. */

  x=y->left;
  y->left=x->right;

  if (nil != x->right)  x->right->parent=y; /*used to use sentinel here */
  /* and do an unconditional assignment instead of testing for nil */

  /* instead of checking if x->parent is the root as in the book, we */
  /* count on the root sentinel to implicitly take care of this case */
  x->parent=y->parent;
  if( y == y->parent->left) {
    y->parent->left=x;
  } else {
    y->parent->right=x;
  }
  x->right=y;
  y->parent=x;

  /*Update rcoords*/
  max_rcoord = ((node_info*)(x->info))->max_node_rcoord;
  if( (x->left != nil) && (((node_info*)x->left->info)->max_subtr_rcoord > max_rcoord))
	max_rcoord = ((node_info*)x->left->info)->max_subtr_rcoord;
  if( (x->right != nil) && (((node_info*)x->right->info)->max_subtr_rcoord > max_rcoord))
	max_rcoord = ((node_info*)x->right->info)->max_subtr_rcoord;
  ((node_info*)x->info)->max_subtr_rcoord = max_rcoord;

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not red in RightRotate");
#endif
}

/***********************************************************************/
/*  FUNCTION:  TreeInsertHelp  */
/**/
/*  INPUTS:  tree is the tree to insert into and z is the node to insert */
/**/
/*  OUTPUT:  none */
/**/
/*  Modifies Input:  tree, z */
/**/
/*  EFFECTS:  Inserts z into the tree as if it were a regular binary tree */
/*            using the algorithm described in _Introduction_To_Algorithms_ */
/*            by Cormen et al.  This funciton is only intended to be called */
/*            by the RBTreeInsert function and not by the user */
/***********************************************************************/

void TreeInsertHelp(rb_red_blk_tree* tree, rb_red_blk_node* z) {
  /*  This function should only be called by InsertRBTree (see above) */
  rb_red_blk_node* x;
  rb_red_blk_node* y;
  rb_red_blk_node* nil=tree->nil;

  z->left=z->right=nil;
  y=tree->root;
  x=tree->root->left;
  while( x != nil) {
    y=x;
    if (1 == tree->Compare(x->key,z->key)) { /* x.key > z.key */
      x=x->left;
    } else { /* x,key <= z.key */
      x=x->right;
    }
  }
  z->parent=y;
  if ( (y == tree->root) ||
       (1 == tree->Compare(y->key,z->key))) { /* y.key > z.key */
    y->left=z;
  } else {
    y->right=z;
  }

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not red in TreeInsertHelp");
#endif
}


/* This function inserts a node in an interval tree representing
 * a multi-dimensional block (array). Because the blocks are multi-dimensional, 
 * we embed trees (one for each dimension).  
 * A tree orders start coordinate of the block 
 * in a given dimension, i.e. each node in the tree has a pointer to the tree that 
 * orders the coordinates in the next dimension. E.g. for a 2D blocks with coorner 
 * coordinate strt, the tree on level-0 will hold values strt[0] for each block. 
 * Each node in tree on level-0 will then have a pointer to the tree on level-1
 * that orders coordinates strt[1] for each block.
 * Each node of the interval tree of a given level is augmented with two values:
 *  maximum right end coordinate of the blocks associated with this node as well as 
 * the maximum right coordinate of this node's subtree. This is later used to find 
 * the block that contains the query subblock.
 * */
void TreeInsertHelpVer2(rb_red_blk_tree* tree, rb_red_blk_node* z, insert_info *info) {
 
  rb_red_blk_node* x, *insert_node;
  MPI_Offset max_rcoord;
  rb_red_blk_node* y;
  rb_red_blk_node* nil=tree->nil;
  int cmp;
  int cur_dim = info->cur_dim;
  int ndims = info->ndims;
  void *key = (void*)(&(info->blck->start[cur_dim]));
  int exist = 0, update_subtr_rcoord = 0;
  
  //assert(cur_dim >= 0 && cur_dim < ndims);
  DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Insert key %lld, cur_dim %d, nnodes %lu\n", *(MPI_Offset*)key, cur_dim, tree->nnodes);
  /*If cur_dim == ndims - 1, we are at the highest dimension, so we insert node z, 
   * otherwise, create an intermediate node for this tree of level cur_dim and 
   * then proceed to the tree in the dimension cur_dim+1. 
   * After each insertion we need to update node colors and augmented info of 
   * the nodes on the path from the inserted node to the root. */
  
  /*Find where to insert the node*/
  y=tree->root;
  x=tree->root->left;
  while( x != nil) {
    y=x;
    cmp = tree->Compare(x->key,key);
    if (cmp == 1) { /* x.key > key */
      x=x->left;
      DTF_DBG(VERBOSE_RB_TREE_LEVEL, "->left");
    } else if (cmp == -1){ /* x.key < key */
      x=x->right;
      DTF_DBG(VERBOSE_RB_TREE_LEVEL, "->left");
    } else {
		DTF_DBG(VERBOSE_RB_TREE_LEVEL, "->exist");
		//DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Trying to insert key %lld to node with key ", *(MPI_Offset*)key);
		//tree->PrintKey(x->key);
		//DTF_DBG(VERBOSE_RB_TREE_LEVEL,"\n");
		
		//Node already exists, update max values and go to next dimension
		exist = 1;
		insert_node = x;
		max_rcoord = info->blck->start[cur_dim] + info->blck->count[cur_dim] - 1;
		
		if(cur_dim == ndims - 1){//blocks overlap (overwriting data)
			int i;
			//allow only exact same data overwrite
			if(max_rcoord != ((node_info*)(x->info))->max_node_rcoord){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: process tries to overwrite data, but the new block only partially overlaps with the previously written block. We don't allow that.");
				MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
			}
			if(info->blck->rank != ((node_info*)(x->info))->blck->rank){
				DTF_DBG(VERBOSE_ALL_LEVEL, "DTF Warning: Process %d overwrites data previously written by process %d", info->blck->rank, ((node_info*)(x->info))->blck->rank);
				((node_info*)(x->info))->blck->rank = info->blck->rank;
			}
			for (i=0; i < info->ndims; i++)
				DTF_DBG(VERBOSE_ALL_LEVEL, "old: %lld->%lld new: %lld->%lld",((node_info*)(x->info))->blck->start[i], ((node_info*)(x->info))->blck->count[i], info->blck->start[i], info->blck->count[i]); 

		} 
		
		if(max_rcoord > ((node_info*)(x->info))->max_node_rcoord ){
			((node_info*)(x->info))->max_node_rcoord = max_rcoord;
			
			if(max_rcoord > ((node_info*)(x->info))->max_subtr_rcoord){
				((node_info*)(x->info))->max_subtr_rcoord = max_rcoord;
				update_subtr_rcoord = 1;
			}
		}
			
		break;
	}
  }
  
  if(!exist){ 
	  if(cur_dim == ndims - 1){
		z->left=z->right=z->parent=nil;
		insert_node = z;
		//DTF_DBG(VERBOSE_DBG_LEVEL,"Last dim, insert node %p\n", insert_node);
	  } else {
		 
		/*new intermediate node for tree in dimension cur_dim*/
		insert_node = (rb_red_blk_node*)SafeMalloc(sizeof(rb_red_blk_node));
		insert_node->key = key;
  
		node_info *ninfo = (node_info*)SafeMalloc(sizeof(node_info));
	    ninfo->max_node_rcoord = info->blck->start[cur_dim]+info->blck->count[cur_dim]-1;
	    ninfo->max_subtr_rcoord = ninfo->max_node_rcoord;
	    /*Nodes in a tree not of the highest dimension represent not a block but 
	     * a list of blocks whose coordinates overlap in this dimension*/
	    ninfo->blck = NULL; 
	    ninfo->next_dim_tree = RBTreeCreateBlocks(tree->Compare, tree->DestroyKey, 
											tree->DestroyInfo, tree->PrintKey, 
											tree->PrintInfo, cur_dim+1);
	    insert_node->info = (void*)ninfo;
	    insert_node->left=insert_node->right=insert_node->parent=nil;
	    DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Created new intermediate node %p\n", (void*)insert_node); 
	  }
	  nnodes++;
	  tree->nnodes++;
	  //DTF_DBG(VERBOSE_DBG_LEVEL,"node left %p right %p parent %p\n",insert_node->left, insert_node->right, insert_node->parent); 
	  //if(tree->root == nil){
		 //DTF_DBG(VERBOSE_DBG_LEVEL,"Insert to the root\n");
		 //DTF_DBG(VERBOSE_DBG_LEVEL,"nil %p, root %p, left %p, right %p, parent %p\n", nil, tree->root, tree->root->left, tree->root->right, tree->root->parent);
		 //insert_node->red = 0;
		 //tree->root = insert_node;
		 //DTF_DBG(VERBOSE_DBG_LEVEL,"nil %p, root %p, left %p, right %p, parent %p\n", nil, tree->root, tree->root->left, tree->root->right, tree->root->parent);
		 //assert(insert_node->left == nil);
		 //assert(insert_node->right == nil);
	  //} else{
	  insert_node->parent=y;
	  if ( (y == tree->root) ||
		   (1 == tree->Compare(y->key,insert_node->key))) { /* y.key > z.key */
		   y->left=insert_node;
	  } else {
		y->right=insert_node;
	  }
	 // }
 
#ifdef DEBUG_ASSERT
	  Assert(!tree->nil->red,"nil not red in TreeInsertHelp");
#endif
	  //if(insert_node != tree->root){
		  //DTF_DBG(VERBOSE_DBG_LEVEL,"not a root %p (root %p)\n", insert_node, tree->root);
		  //Balance the tree
	  DTF_DBG(VERBOSE_RB_TREE_LEVEL, "Begin balance");
	  x = insert_node;
	  x->red=1;
	  while(x->parent->red) { /* use sentinel instead of checking for root */
		if (x->parent == x->parent->parent->left) {
		  y=x->parent->parent->right;
		  if (y->red) {
			x->parent->red=0;
			y->red=0;
			x->parent->parent->red=1;
			x=x->parent->parent;
		  } else {
			if (x == x->parent->right) {
			  x=x->parent;
			  LeftRotate(tree,x);
			}
			x->parent->red=0;
			x->parent->parent->red=1;
			RightRotate(tree,x->parent->parent);
		  }
		} else { /* case for x->parent == x->parent->parent->right */
		  y=x->parent->parent->left;
		  if (y->red) {
			x->parent->red=0;
			y->red=0;
			x->parent->parent->red=1;
			x=x->parent->parent;
		  } else {
			if (x == x->parent->left) {
			  x=x->parent;
			  RightRotate(tree,x);
			}
			x->parent->red=0;
			x->parent->parent->red=1;
			LeftRotate(tree,x->parent->parent);
		  }
		}
	  }
	  tree->root->left->red=0;
		DTF_DBG(VERBOSE_RB_TREE_LEVEL, "End balance");
	}

  if(!exist || update_subtr_rcoord){
	  //update max rcoords on the path back to the parent
	  x = insert_node;
	  while(x != tree->root){
		  assert(x->info != NULL);
		  max_rcoord = ((node_info*)(x->info))->max_node_rcoord;
		  if( (x->left != nil) && (((node_info*)x->left->info)->max_subtr_rcoord > max_rcoord))
			max_rcoord = ((node_info*)x->left->info)->max_subtr_rcoord;
		  if( (x->right != nil) && (((node_info*)x->right->info)->max_subtr_rcoord > max_rcoord))
			max_rcoord = ((node_info*)x->right->info)->max_subtr_rcoord;
		  ((node_info*)x->info)->max_subtr_rcoord = max_rcoord;
		  x = x->parent;
	  }
	  DTF_DBG(VERBOSE_RB_TREE_LEVEL, "Updated max_subtr_rcoord");
   }
  
  
  if(cur_dim != ndims - 1){
	//goto higher dimension tree
	info->cur_dim++;
  	TreeInsertHelpVer2( ((node_info*)(insert_node->info))->next_dim_tree, z, info );
  }

}
/*  Before calling Insert RBTree the node x should have its key set */

/***********************************************************************/
/*  FUNCTION:  RBTreeInsert */
/**/
/*  INPUTS:  tree is the red-black tree to insert a node which has a key */
/*           pointed to by key and info pointed to by info.  */
/**/
/*  OUTPUT:  This function returns a pointer to the newly inserted node */
/*           which is guarunteed to be valid until this node is deleted. */
/*           What this means is if another data structure stores this */
/*           pointer then the tree does not need to be searched when this */
/*           is to be deleted. */
/**/
/*  Modifies Input: tree */
/**/
/*  EFFECTS:  Creates a node node which contains the appropriate key and */
/*            info pointers and inserts it into the tree. */
/***********************************************************************/

rb_red_blk_node * RBTreeInsert(rb_red_blk_tree* tree, void* key, void* info) {
  rb_red_blk_node * y;
  rb_red_blk_node * x;
  rb_red_blk_node * newNode;
  
  if(tree->root == tree->nil)
	assert(0); //modified create tree function so that root is nil when tree created 

  x=(rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  x->key=key;
  x->info=info;

  TreeInsertHelp(tree,x);
  newNode=x;
  x->red=1;
  while(x->parent->red) { /* use sentinel instead of checking for root */
    if (x->parent == x->parent->parent->left) {
      y=x->parent->parent->right;
      if (y->red) {
	x->parent->red=0;
	y->red=0;
	x->parent->parent->red=1;
	x=x->parent->parent;
      } else {
	if (x == x->parent->right) {
	  x=x->parent;
	  LeftRotate(tree,x);
	}
	x->parent->red=0;
	x->parent->parent->red=1;
	RightRotate(tree,x->parent->parent);
      }
    } else { /* case for x->parent == x->parent->parent->right */
      y=x->parent->parent->left;
      if (y->red) {
	x->parent->red=0;
	y->red=0;
	x->parent->parent->red=1;
	x=x->parent->parent;
      } else {
	if (x == x->parent->left) {
	  x=x->parent;
	  RightRotate(tree,x);
	}
	x->parent->red=0;
	x->parent->parent->red=1;
	LeftRotate(tree,x->parent->parent);
      }
    }
  }
  tree->root->left->red=0;
  return(newNode);

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not red in RBTreeInsert");
  Assert(!tree->root->red,"root not red in RBTreeInsert");
#endif
}

rb_red_blk_node * RBTreeInsertBlock(rb_red_blk_tree* tree, void* ins_info) 
{
  rb_red_blk_node * x;
  rb_red_blk_node * newNode;
  insert_info *info = (insert_info*)ins_info;
  node_info *ninfo;
  int ndims = info->ndims;
  
  /*this block will be inserted in the tree on the 
   * level of the last dimension, i.e., this is a leaf node*/
  x=(rb_red_blk_node*) SafeMalloc(sizeof(rb_red_blk_node));
  nnodes++;
  x->key=(void*)(&(info->blck->start[info->ndims-1]));
  
  ninfo = (node_info*)SafeMalloc(sizeof(node_info));
  /*Because this will be a leaf node, max value of the 
   * right coordinate is the max value for this block*/
  ninfo->max_node_rcoord = info->blck->start[ndims-1]+info->blck->count[ndims-1]-1;
  ninfo->max_subtr_rcoord = ninfo->max_node_rcoord;
  ninfo->blck = info->blck;
  ninfo->next_dim_tree = NULL;
  x->info = ninfo;
  
  //x->left=x->right=x->parent=nil; //this init should be done inside the respective tree

  TreeInsertHelpVer2(tree,x, info);
  newNode=x;
  
#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not red in RBTreeInsert");
  Assert(!tree->root->red,"root not red in RBTreeInsert");
#endif
return(newNode);
}

/***********************************************************************/
/*  FUNCTION:  TreeSuccessor  */
/**/
/*    INPUTS:  tree is the tree in question, and x is the node we want the */
/*             the successor of. */
/**/
/*    OUTPUT:  This function returns the successor of x or NULL if no */
/*             successor exists. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
/***********************************************************************/

rb_red_blk_node* TreeSuccessor(rb_red_blk_tree* tree,rb_red_blk_node* x) {
  rb_red_blk_node* y;
  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* root=tree->root;

  if (nil != (y = x->right)) { /* assignment to y is intentional */
    while(y->left != nil) { /* returns the minium of the right subtree of x */
      y=y->left;
    }
    return(y);
  } else {
    y=x->parent;
    while(x == y->right) { /* sentinel used instead of checking for nil */
      x=y;
      y=y->parent;
    }
    if (y == root) return(nil);
    return(y);
  }
}

/***********************************************************************/
/*  FUNCTION:  Treepredecessor  */
/**/
/*    INPUTS:  tree is the tree in question, and x is the node we want the */
/*             the predecessor of. */
/**/
/*    OUTPUT:  This function returns the predecessor of x or NULL if no */
/*             predecessor exists. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
/***********************************************************************/

rb_red_blk_node* TreePredecessor(rb_red_blk_tree* tree, rb_red_blk_node* x) {
  rb_red_blk_node* y;
  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* root=tree->root;

  if (nil != (y = x->left)) { /* assignment to y is intentional */
    while(y->right != nil) { /* returns the maximum of the left subtree of x */
      y=y->right;
    }
    return(y);
  } else {
    y=x->parent;
    while(x == y->left) {
      if (y == root) return(nil);
      x=y;
      y=y->parent;
    }
    return(y);
  }
}

/***********************************************************************/
/*  FUNCTION:  InorderTreePrint */
/**/
/*    INPUTS:  tree is the tree to print and x is the current inorder node */
/**/
/*    OUTPUT:  none  */
/**/
/*    EFFECTS:  This function recursively prints the nodes of the tree */
/*              inorder using the PrintKey and PrintInfo functions. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:    This function should only be called from RBTreePrint */
/***********************************************************************/

void InorderTreePrint(rb_red_blk_tree* tree, rb_red_blk_node* x) {

  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* root=tree->root;
  if (x != tree->nil) {
    InorderTreePrint(tree,x->left);
    DTF_DBG(VERBOSE_DBG_LEVEL,"info=");
    tree->PrintInfo(x->info);
    DTF_DBG(VERBOSE_DBG_LEVEL,"  key=");
    tree->PrintKey(x->key);
    DTF_DBG(VERBOSE_DBG_LEVEL,"  l->key=");
    if( x->left == nil) DTF_DBG(VERBOSE_DBG_LEVEL,"NULL"); else tree->PrintKey(x->left->key);
    DTF_DBG(VERBOSE_DBG_LEVEL,"  r->key=");
    if( x->right == nil) DTF_DBG(VERBOSE_DBG_LEVEL,"NULL"); else tree->PrintKey(x->right->key);
    DTF_DBG(VERBOSE_DBG_LEVEL,"  p->key=");
    if( x->parent == root) DTF_DBG(VERBOSE_DBG_LEVEL,"NULL"); else tree->PrintKey(x->parent->key);
    DTF_DBG(VERBOSE_DBG_LEVEL,"  red=%i\n",x->red);
    InorderTreePrint(tree,x->right);
  }
}

/*Prints an embeded interval tree*/
void InorderTreePrintVer2(rb_red_blk_tree* tree, rb_red_blk_node* x) {
  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* root=tree->root;
  
  if (x != tree->nil) {
    InorderTreePrintVer2(tree,x->left);

    printf("key=");
    tree->PrintKey(x->key);
    printf("  l->key=");
    
    if( x->left == nil) printf("nil"); else tree->PrintKey(x->left->key);
    printf("  r->key=");
    if( x->right == nil) printf("nil"); else tree->PrintKey(x->right->key);
    printf("  p->key=");
    if( x->parent == root) printf("nil"); else tree->PrintKey(x->parent->key);
    //DTF_DBG(VERBOSE_DBG_LEVEL,"  red=%i\n",x->red);
    printf("  info=");
    tree->PrintInfo(x->info);
    
    InorderTreePrintVer2(tree,x->right);
  } else 
	DTF_DBG(VERBOSE_DBG_LEVEL,"nil");
}


/***********************************************************************/
/*  FUNCTION:  TreeDestHelper */
/**/
/*    INPUTS:  tree is the tree to destroy and x is the current node */
/**/
/*    OUTPUT:  none  */
/**/
/*    EFFECTS:  This function recursively destroys the nodes of the tree */
/*              postorder using the DestroyKey and DestroyInfo functions. */
/**/
/*    Modifies Input: tree, x */
/**/
/*    Note:    This function should only be called by RBTreeDestroy */
/***********************************************************************/

void TreeDestHelper(rb_red_blk_tree* tree, rb_red_blk_node* x) {
  rb_red_blk_node* nil=tree->nil;
  if (x != nil) {
    TreeDestHelper(tree,x->left);
    TreeDestHelper(tree,x->right);
    tree->DestroyKey(x->key);
    tree->DestroyInfo(x->info);
    free(x);
    assert(tree->nnodes > 0);
    tree->nnodes--;
  }
}


/***********************************************************************/
/*  FUNCTION:  RBTreeDestroy */
/**/
/*    INPUTS:  tree is the tree to destroy */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  Destroys the key and frees memory */
/**/
/*    Modifies Input: tree */
/**/
/***********************************************************************/

void RBTreeDestroy(rb_red_blk_tree* tree) {
  //printf("nnodes %lu\n", tree->nnodes);
  TreeDestHelper(tree,tree->root->left);
  free(tree->root);
  tree->nnodes--;
  free(tree->nil);
  tree->nnodes--;
  
  assert(tree->nnodes == 0);
  free(tree);
}


/***********************************************************************/
/*  FUNCTION:  RBTreePrint */
/**/
/*    INPUTS:  tree is the tree to print */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  This function recursively prints the nodes of the tree */
/*             inorder using the PrintKey and PrintInfo functions. */
/**/
/*    Modifies Input: none */
/**/
/***********************************************************************/

void RBTreePrint(rb_red_blk_tree* tree) {
  InorderTreePrint(tree,tree->root->left);
}

void RBTreePrintBlocks(rb_red_blk_tree *tree) {
	InorderTreePrintVer2(tree,tree->root->left);
}


/***********************************************************************/
/*  FUNCTION:  RBExactQuery */
/**/
/*    INPUTS:  tree is the tree to print and q is a pointer to the key */
/*             we are searching for */
/**/
/*    OUTPUT:  returns the a node with key equal to q.  If there are */
/*             multiple nodes with key equal to q this function returns */
/*             the one highest in the tree */
/**/
/*    Modifies Input: none */
/**/
/***********************************************************************/

rb_red_blk_node* RBExactQuery(rb_red_blk_tree* tree, void* q) {
  rb_red_blk_node* x=tree->root->left;
  rb_red_blk_node* nil=tree->nil;
  int compVal;
  if (x == nil) return(0);
  compVal=tree->Compare(x->key,(int*) q);
  while(0 != compVal) {/*assignemnt*/
    if (1 == compVal) { /* x->key > q */
      x=x->left;
    } else {
      x=x->right;
    }
    if ( x == nil) return(0);
    compVal=tree->Compare(x->key,(int*) q);
  }
  return(x);
}

rb_red_blk_node* RBExactQueryPoint(rb_red_blk_tree* tree, MPI_Offset* start, int ndims) {
  rb_red_blk_node* x=tree->root->left;
  rb_red_blk_node* nil=tree->nil;
  if (x == nil) return(0);
  return FindBlockPoint(tree, x, start, ndims);
  
}

rb_red_blk_node* RBExactQueryOverlap(rb_red_blk_tree* tree, MPI_Offset* start, MPI_Offset* count, int ndims) {
  rb_red_blk_node* x=tree->root->left;
  rb_red_blk_node* nil=tree->nil;
  if (x == nil) return(0);
  return FindBlockOverlap(tree, x, start, count, ndims);
}

/*find block that includes the [start] coordinate*/
rb_red_blk_node* FindBlockPoint(rb_red_blk_tree* tree, rb_red_blk_node *subtr, MPI_Offset* start, int ndims)
{
	node_info *info;
	rb_red_blk_node *node = NULL, *nil = tree->nil;
	int cur_dim = tree->cur_dim;
	MPI_Offset coord = start[cur_dim];
	
    if(subtr == nil) return 0;
    
    info = (node_info*)subtr->info;
    DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Looking for %lld\n", coord);
    
    DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Cur node %lld, max rcoord %lld, max tree rcoord %lld\n", 
			*(MPI_Offset*)subtr->key, info->max_node_rcoord, info->max_subtr_rcoord);
    
    if( coord >= *(MPI_Offset*)(subtr->key) && 
        coord <= info->max_node_rcoord ){
			
			if(coord <= info->max_node_rcoord){
			    if(cur_dim == ndims - 1){
					DTF_DBG(VERBOSE_RB_TREE_LEVEL,"found\n");
					//found node
					return subtr;
				} else {
					//DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go deeper\n");
					//go to next dim
					node = FindBlockPoint( info->next_dim_tree, info->next_dim_tree->root->left /*true root*/, start, ndims);
				}
			}
			if(node != NULL) return node;
	}		
	
	if( subtr->left != nil && 
		//coord >= *(MPI_Offset*)(subtr->left->key) &&
		coord <= ((node_info*)subtr->left->info)->max_subtr_rcoord ){
		DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go left\n");
		//investigate left subtree
		node = FindBlockPoint( tree, subtr->left, start, ndims);
		if(node != NULL) return node;
	} 
	
	if( subtr->right != nil && 
		//coord >= *(MPI_Offset*)(subtr->right->key) &&
		coord <= ((node_info*)subtr->right->info)->max_subtr_rcoord ){
		DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go right\n");
		//investigate right subtree
		node = FindBlockPoint( tree, subtr->right, start, ndims);  
		if(node != NULL) return node;
	} 
		
	DTF_DBG(VERBOSE_RB_TREE_LEVEL,"nope\n");
	return node;
}

/*Find block that overlaps in any way with quired block*/
rb_red_blk_node* FindBlockOverlap(rb_red_blk_tree* tree, rb_red_blk_node *subtr, MPI_Offset* start, MPI_Offset* count, int ndims)
{
	node_info *info;
	rb_red_blk_node *node = NULL, *nil = tree->nil;
	int cur_dim = tree->cur_dim;
	MPI_Offset lcoord = start[cur_dim];
	MPI_Offset rcoord = start[cur_dim] + count[cur_dim] - 1;
    if(subtr == nil) return 0;
    
    info = (node_info*)subtr->info;
    DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Looking for [%lld, %lld]\n", lcoord, rcoord);
    
    DTF_DBG(VERBOSE_RB_TREE_LEVEL,"Cur node %lld, max rcoord %lld, max tree rcoord %lld\n", 
			*(MPI_Offset*)subtr->key, info->max_node_rcoord, info->max_subtr_rcoord);
    
    if( (lcoord >= *(MPI_Offset*)(subtr->key) && lcoord <= info->max_node_rcoord) ||
        (lcoord < *(MPI_Offset*)(subtr->key) &&  rcoord >= *(MPI_Offset*)(subtr->key))){
			
			
			if(cur_dim == ndims - 1){
				DTF_DBG(VERBOSE_RB_TREE_LEVEL,"found\n");
				//found node
				return subtr;
			} else {
				//DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go deeper\n");
				//go to next dim
				node = FindBlockOverlap( info->next_dim_tree, info->next_dim_tree->root->left /*true root*/, start, count, ndims);
			}
			
			if(node != NULL) return node;
	}		
	
	if( subtr->left != nil && 
		//coord >= *(MPI_Offset*)(subtr->left->key) &&
		lcoord <= ((node_info*)subtr->left->info)->max_subtr_rcoord){
		DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go left\n");
		//investigate left subtree
		node = FindBlockOverlap( tree, subtr->left, start,count, ndims);
		if(node != NULL) return node;
	} 
	
	if( subtr->right != nil && 
		rcoord >= *(MPI_Offset*)(subtr->right->key) &&
		lcoord <= ((node_info*)subtr->right->info)->max_subtr_rcoord ){
		DTF_DBG(VERBOSE_RB_TREE_LEVEL,"go right\n");
		//investigate right subtree
		node = FindBlockOverlap( tree, subtr->right, start, count, ndims);  
		if(node != NULL) return node;
	} 
		
	DTF_DBG(VERBOSE_RB_TREE_LEVEL,"nope\n");
	return node;
}

/***********************************************************************/
/*  FUNCTION:  RBDeleteFixUp */
/**/
/*    INPUTS:  tree is the tree to fix and x is the child of the spliced */
/*             out node in RBTreeDelete. */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  Performs rotations and changes colors to restore red-black */
/*             properties after a node is deleted */
/**/
/*    Modifies Input: tree, x */
/**/
/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
/***********************************************************************/

void RBDeleteFixUp(rb_red_blk_tree* tree, rb_red_blk_node* x) {
  rb_red_blk_node* root=tree->root->left;
  rb_red_blk_node* w;

  while( (!x->red) && (root != x)) {
    if (x == x->parent->left) {
      w=x->parent->right;
      if (w->red) {
	w->red=0;
	x->parent->red=1;
	LeftRotate(tree,x->parent);
	w=x->parent->right;
      }
      if ( (!w->right->red) && (!w->left->red) ) {
	w->red=1;
	x=x->parent;
      } else {
	if (!w->right->red) {
	  w->left->red=0;
	  w->red=1;
	  RightRotate(tree,w);
	  w=x->parent->right;
	}
	w->red=x->parent->red;
	x->parent->red=0;
	w->right->red=0;
	LeftRotate(tree,x->parent);
	x=root; /* this is to exit while loop */
      }
    } else { /* the code below is has left and right switched from above */
      w=x->parent->left;
      if (w->red) {
	w->red=0;
	x->parent->red=1;
	RightRotate(tree,x->parent);
	w=x->parent->left;
      }
      if ( (!w->right->red) && (!w->left->red) ) {
	w->red=1;
	x=x->parent;
      } else {
	if (!w->left->red) {
	  w->right->red=0;
	  w->red=1;
	  LeftRotate(tree,w);
	  w=x->parent->left;
	}
	w->red=x->parent->red;
	x->parent->red=0;
	w->left->red=0;
	RightRotate(tree,x->parent);
	x=root; /* this is to exit while loop */
      }
    }
  }
  x->red=0;

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not black in RBDeleteFixUp");
#endif
}


/***********************************************************************/
/*  FUNCTION:  RBDelete */
/**/
/*    INPUTS:  tree is the tree to delete node z from */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  Deletes z from tree and frees the key and info of z */
/*             using DestoryKey and DestoryInfo.  Then calls */
/*             RBDeleteFixUp to restore red-black properties */
/**/
/*    Modifies Input: tree, z */
/**/
/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
/***********************************************************************/

void RBDelete(rb_red_blk_tree* tree, rb_red_blk_node* z){
  rb_red_blk_node* y;
  rb_red_blk_node* x;
  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* root=tree->root;

  y= ((z->left == nil) || (z->right == nil)) ? z : TreeSuccessor(tree,z);
  x= (y->left == nil) ? y->right : y->left;
  if (root == (x->parent = y->parent)) { /* assignment of y->p to x->p is intentional */
    root->left=x;
  } else {
    if (y == y->parent->left) {
      y->parent->left=x;
    } else {
      y->parent->right=x;
    }
  }
  if (y != z) { /* y should not be nil in this case */

#ifdef DEBUG_ASSERT
    Assert( (y!=tree->nil),"y is nil in RBDelete\n");
#endif
    /* y is the node to splice out and x is its child */

    if (!(y->red)) RBDeleteFixUp(tree,x);

    tree->DestroyKey(z->key);
    tree->DestroyInfo(z->info);
    y->left=z->left;
    y->right=z->right;
    y->parent=z->parent;
    y->red=z->red;
    z->left->parent=z->right->parent=y;
    if (z == z->parent->left) {
      z->parent->left=y;
    } else {
      z->parent->right=y;
    }
    free(z);
  } else {
    tree->DestroyKey(y->key);
    tree->DestroyInfo(y->info);
    if (!(y->red)) RBDeleteFixUp(tree,x);
    free(y);
  }

#ifdef DEBUG_ASSERT
  Assert(!tree->nil->red,"nil not black in RBDelete");
#endif
}


/***********************************************************************/
/*  FUNCTION:  RBDEnumerate */
/**/
/*    INPUTS:  tree is the tree to look for keys >= low */
/*             and <= high with respect to the Compare function */
/**/
/*    OUTPUT:  stack containing pointers to the nodes between [low,high] */
/**/
/*    Modifies Input: none */
/***********************************************************************/

stk_stack* RBEnumerate(rb_red_blk_tree* tree, void* low, void* high) {
  stk_stack* enumResultStack;
  rb_red_blk_node* nil=tree->nil;
  rb_red_blk_node* x=tree->root->left;
  rb_red_blk_node* lastBest=nil;

  enumResultStack=StackCreate();
  while(nil != x) {
    if ( 1 == (tree->Compare(x->key,high)) ) { /* x->key > high */
      x=x->left;
    } else {
      lastBest=x;
      x=x->right;
    }
  }
  while ( (lastBest != nil) && (1 != tree->Compare(low,lastBest->key))) {
    StackPush(enumResultStack,lastBest);
    lastBest=TreePredecessor(tree,lastBest);
  }
  return(enumResultStack);
}

void rb_print_stats()
{
	if(ntrees > 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL,"Allocd %llu trees (%llu bytes), %llu nodes (%llu bytes)\n", ntrees, ntrees*sizeof(rb_red_blk_tree),  
		        nnodes, nnodes*sizeof(rb_red_blk_node)*sizeof(node_info));
}


/***************************FUNCTIONS USED IN DTF*******************************/

void bl_print(block_t *_bl, int ndims)
{
  //printf(" node %p\n", (void*)_bl);
  int i;
  printf("([");
  for(i = 0; i < ndims; i++)
	printf("%lld,", _bl->start[i]);
  printf("],["); 
  for(i = 0; i < ndims; i++)
	printf("%lld,", _bl->count[i]);
  printf("]) \n");
}

/*API for handling rb_tree in write_db_item*/
void rb_destroy_node_info(void* ninfo)
{
	node_info *info = (node_info*)ninfo;
	if(info->next_dim_tree != NULL)
		RBTreeDestroy(info->next_dim_tree);
	
	if(info->blck != NULL){
		free(info->blck->start);
		free(info->blck->count);
		free(info->blck);
	}
	
	free(info);
}

int rb_key_cmp(const void *a, const void *b)
{
 
	if(*(MPI_Offset*)a > *(MPI_Offset*)b)
		return 1;
	else if (*(MPI_Offset*)a < *(MPI_Offset*)b)
		return -1;
  return 0;
}

void rb_print_key(const void *key)
{
	printf("%lld", *(MPI_Offset*)key);
}

void rb_print_info(void *ninfo)
{
	node_info *info = (node_info*)ninfo;
	printf(" max rcoord: %lld", info->max_node_rcoord);
	printf(" max subtr rcoord: %lld", info->max_subtr_rcoord);
	
	if(info->blck != NULL){
		bl_print(info->blck, 2);
	} else 
		printf("\n");
	
	if(info->next_dim_tree != NULL){
		DTF_DBG(VERBOSE_DBG_LEVEL,"Tree for dim %d---->\n", info->next_dim_tree->cur_dim);
		RBTreePrintBlocks(info->next_dim_tree);
		DTF_DBG(VERBOSE_DBG_LEVEL,"<-----Tree for dim %d\n", info->next_dim_tree->cur_dim);
	}
}

block_t *rb_find_interval_point(rb_red_blk_tree* tree, MPI_Offset* start, int ndims)
{
	rb_red_blk_node *find=RBExactQueryPoint(tree, start, ndims);
	
	if(find != NULL)
		return ((node_info*)find->info)->blck;
	else return NULL;
}

block_t *rb_find_interval_overlap(rb_red_blk_tree* tree, MPI_Offset* start, MPI_Offset* count, int ndims)
{
	rb_red_blk_node *find=RBExactQueryOverlap(tree, start, count, ndims);
	
	if(find != NULL)
		return ((node_info*)find->info)->blck;
	else return NULL;
}

block_t *rb_find_block(rb_red_blk_tree* tree, MPI_Offset* start, MPI_Offset* count, int ndims)
{
	return rb_find_interval_point(tree, start, ndims);
	//return rb_find_interval_overlap(tree, start, count, ndims);
}

