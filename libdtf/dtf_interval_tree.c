#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dtf_interval_tree.h"
#include "dtf_util.h"

static size_t memuse = 0;
static unsigned long nnodes = 0, ntrees = 0;

static int cmp_intervals(MPI_Offset la, MPI_Offset ra, MPI_Offset lb, MPI_Offset rb)
{
	if(lb >= ra) return 1; //interval b lies completely to the right from a
	if(rb <= la) return -1; //interval b lies completely to the left from a
	return 0; //they overlap
}

static interval_node_t* create_node(MPI_Offset lkey, MPI_Offset rkey, block_t *bl, int cur_dim, int ndims, int node_id)
{
	interval_node_t *node = dtf_malloc(sizeof(interval_node_t));
	memuse += sizeof(interval_node_t);
	node->id = node_id;
	node_id++;
	node->lkey = lkey;
	node->rkey = rkey;
	node->bl_refs = NULL;
	node->height = 1;
	node->left=node->right=node->parent=NULL;
	node->next_dim_tree = NULL;
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Created node %d for interval [%lld, %lld) in dim %d", node->id, lkey, rkey, cur_dim);
	nnodes++;
	
	return node;
}

static void delete_node(interval_node_t *node)
{
	block_ref_t* ref;
	
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"delete node %d", node->id);
	
	if(node->next_dim_tree != NULL){
		IntervalTreeDestroy(node->next_dim_tree);
		node->next_dim_tree = NULL;
	}
	
	//~ if(node->parent != NULL){
		//~ if(node == node->parent->left) 
			//~ node->parent->left = NULL;
		//~ else 
			//~ node->parent->right = NULL;
	//~ }
	ref = node->bl_refs;
	while(ref != NULL){
		node->bl_refs = node->bl_refs->next;
		dtf_free(ref, sizeof(block_ref_t));
		ref = node->bl_refs;
	}
	dtf_free(node, sizeof(interval_node_t));
}

static void destroy_subtree(interval_node_t* node)
{
	if(node == NULL) return;
	
	if(node->left != NULL)
		destroy_subtree(node->left);
	if(node->right != NULL)
		destroy_subtree(node->right);
	
	delete_node(node);
}

static int inline max(int a, int b)
{
	return a>b?a:b;
}

static int inline height(interval_node_t* node)
{
	return (node == NULL) ? 0: node->height;
}

static int inline balance(interval_node_t* node)
{
	return (node == NULL) ? 0: height(node->left) - height(node->right);
}

static interval_node_t* rotate_right(interval_node_t* pivot)
{
	interval_node_t *ltree = pivot->left;
	
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"rotate right at node %d", pivot->id);
	
	pivot->left = ltree->right;
	ltree->right = pivot; 
	
	pivot->height = 1 + max(height(pivot->left), height(pivot->right));
	ltree->height = 1 + max(height(ltree->left), height(ltree->right));
	
	return ltree;
}

static interval_node_t *rotate_left(interval_node_t* pivot)
{
	interval_node_t *rtree = pivot->right;
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"rotate left at node %d", pivot->id);
	pivot->right = rtree->left;
	rtree->left = pivot;
	pivot->height = 1 + max(height(pivot->left), height(pivot->right));
	rtree->height = 1 + max(height(rtree->left), height(rtree->right));
	return rtree;
}

static interval_node_t* insert(interval_node_t* node, interval_tree_t* tree,MPI_Offset lkey, MPI_Offset rkey, block_t* bl, int ndims)
{
	int res;
	
	if(node == NULL){
		node = create_node(lkey, rkey, bl, tree->cur_dim, ndims, tree->nnodes);
		tree->nnodes++;
		if(ndims == 0 || tree->cur_dim == ndims - 1){ 
			block_ref_t *ref = dtf_malloc(sizeof(block_ref_t));
			memuse += sizeof(block_ref_t);
			ref->bl = bl;
			ref->next = node->bl_refs;
			node->bl_refs = ref;
			DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Added ref to block %d in node %d", ref->bl->id, node->id);
		} else {
			int dim = tree->cur_dim + 1;
			if(node->next_dim_tree == NULL){
				node->next_dim_tree = IntervalTreeCreate();
				node->next_dim_tree->cur_dim = dim;
			}
			node->next_dim_tree->root = insert(node->next_dim_tree->root, node->next_dim_tree, bl->start[dim], bl->start[dim]+bl->count[dim], bl, ndims);
		}
		return node;
	}
	//DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"trying to insert interval [%lld, %lld) on level %d", lkey, rkey, tree->cur_dim);
	
	res = cmp_intervals(node->lkey, node->rkey, lkey, rkey);
	if(res < 0) {
		DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"go left");
		node->left = insert(node->left, tree, lkey, rkey, bl, ndims);
	} else if( res > 0) {
		DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"go right");
		node->right = insert(node->right, tree, lkey, rkey, bl, ndims);
	} else {
		DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"they overlap");
		//insert nodes for non-overlapping part if any
		if(lkey < node->lkey){
			DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"go left for [%lld, %lld)", lkey, node->lkey);
			node->left = insert(node->left, tree, lkey, node->lkey, bl, ndims);
		}if(rkey > node->rkey){
			DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"go right for [%lld, %lld)", node->rkey, rkey);
			node->right = insert(node->right, tree, node->rkey, rkey, bl, ndims);
		}
		if(tree->cur_dim == ndims - 1){ 
			block_ref_t* ref = dtf_malloc(sizeof(block_ref_t));
			memuse += sizeof(block_ref_t);
			ref->bl = bl;
			assert(node->bl_refs != NULL);
			ref->next = node->bl_refs;
			node->bl_refs = ref;
			DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Added ref to block %d in node %d", ref->bl->id, node->id);
		} else {
			int dim = tree->cur_dim + 1;
			if(node->next_dim_tree == NULL){
				node->next_dim_tree = IntervalTreeCreate();
				node->next_dim_tree->cur_dim = dim;
			}
			node->next_dim_tree->root = insert(node->next_dim_tree->root, node->next_dim_tree, bl->start[dim], bl->start[dim]+bl->count[dim], bl, ndims);
		}
		return node;
	}		
	
	node->height = 1 + max(height(node->left), height(node->right));
	
	/*left left*/
	if(balance(node) > 1 && balance(node->left) >= 1)
		node = rotate_right(node);
	/*left right*/
	if(balance(node) > 1 && balance(node->left) <= -1){
		node->left = rotate_left(node->left);
		node = rotate_right(node);
	}
	/*right right*/
	if(balance(node) < -1 && balance(node->right) <= -1)
		node = rotate_left(node);
	/*right left*/
	if(balance(node) < -1 && balance(node->right) >= 1){
		node->right = rotate_right(node->right);
		node = rotate_left(node);
	}
	
	return node;
		
}


static void print_node(interval_node_t *node)
{
	block_ref_t *ref = node->bl_refs;
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Node %d, interval [%lld, %lld), height %d ", node->id, node->lkey, node->rkey, height(node));
	if(node->parent != NULL) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"parent %d, ", node->parent->id);
	if(node->left != NULL) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"left %d, ", node->left->id);
	if(node->right != NULL) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"right %d, ", node->right->id);
	if(ref != NULL) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"refs: ");
	while(ref != NULL){
		DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"%d, ", ref->bl->id);
		ref = ref->next;
	}
	
	if(node->next_dim_tree != NULL) IntervalTreePrint(node->next_dim_tree);
}

static void print_subtree(interval_node_t *node)
{
	if(node == NULL) return;
	print_node(node);
	if(node->left != NULL) print_subtree(node->left);
	if(node->right != NULL) print_subtree(node->right);
}

static interval_node_t* find_overlap(interval_node_t* node, MPI_Offset* start, MPI_Offset* count, int cur_dim, int ndims)
{
	MPI_Offset lkey, rkey;
	interval_node_t* ret_node = NULL;
	if(node == NULL) return NULL;
	if(ndims == 0) return node;
	
	lkey = start[cur_dim];
	rkey = start[cur_dim] + count[cur_dim];
	
	
	if(lkey < node->lkey){
		ret_node = find_overlap(node->left, start, count, cur_dim, ndims);
	}
	
	//DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "investigate node %d in dim %d\n", node->id, cur_dim);
	
	if( cmp_intervals(node->lkey, node->rkey, lkey, rkey) == 0){
		if(cur_dim == ndims - 1) 
			ret_node = node;
		 else {
			ret_node = find_overlap(node->next_dim_tree->root, start, count, cur_dim+1, ndims);
		}
	} 
	
	if(ret_node != NULL) return ret_node;
	
	if(rkey > node->rkey){
		ret_node = find_overlap(node->right, start, count, cur_dim, ndims);
	}
	
	
	return ret_node;
}

/*-----------------------------------------------------*/

interval_tree_t* IntervalTreeCreate()
{
	interval_tree_t* tree = dtf_malloc(sizeof(interval_tree_t));
	memuse += sizeof(interval_tree_t);
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Create tree");
	tree->root = NULL;
	tree->cur_dim = 0;
	tree->nnodes = 0;
	ntrees++;
	return tree;
}

void IntervalTreeDestroy(interval_tree_t* tree)
{
	if(tree != NULL)
		DTF_DBG(VERBOSE_DBG_LEVEL, "Destroy tree level %d, height %d, nnodes %lu", tree->cur_dim, height(tree->root), tree->nnodes);
	destroy_subtree(tree->root);
	dtf_free(tree, sizeof(interval_tree_t));	
}

void IntervalTreeInsert(interval_tree_t* tree, insert_info_t* info)
{
	MPI_Offset lkey, rkey;
	int i;
	
	if(info->ndims == 0)
		lkey = rkey = 0;
	else {
		lkey = info->bl->start[0];
		rkey = info->bl->start[0] + info->bl->count[0];
	}

	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Insert block (");
	for(i = 0; i < info->ndims; i++) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"%lld ", info->bl->start[i]);
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"),(");
	for(i = 0; i < info->ndims; i++) DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"%lld ", info->bl->count[i]);
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,")");
	
	tree->root = insert(tree->root, tree, lkey, rkey, info->bl, info->ndims);
}

void IntervalTreePrint(interval_tree_t *tree)
{
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL,"Printing tree on level %d", tree->cur_dim);
	print_subtree(tree->root);
}

/*Returns the first block in the tree with which current lock overlaps*/

block_t* IntervalTreeFindOverlap(interval_tree_t* tree, MPI_Offset* start, MPI_Offset* count, int ndims)
{
	int i;
	block_t *bl = NULL;
	
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "Search for block: (");
	for(i=0;i<ndims;i++)DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "%lld,", start[i]);
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "),(");
	for(i=0;i<ndims;i++)DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "%lld,", count[i]);
	DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, ")");
	
	if(ndims == 0) //scalar var
		return tree->root->bl_refs->bl;
	
	interval_node_t* node = find_overlap(tree->root, start, count, 0, ndims);
	if(node == NULL){
		DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "didn't find overlap\n");
		return bl;
	}
	block_ref_t* ref = node->bl_refs;
	while(ref != NULL){
		for(i = 0; i < ndims; i++)
			if(cmp_intervals(ref->bl->start[i], ref->bl->start[i] + ref->bl->count[i], start[i], start[i] + count[i]) != 0) break;

		if(i == ndims){
			bl = ref->bl;
			DTF_DBG(VERBOSE_INTERVAL_TREE_LEVEL, "found overlap with block %d\n", bl->id);
			break;
		}
		ref = ref->next;
	}

	return bl;
}

void IntervalTreePrintMem()
{
	if(memuse > 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL, "Memory use for trees %ld, %lu nodes, %lu trees", memuse, nnodes, ntrees);
}
