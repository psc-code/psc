///\file bintree.h
///This provides an implementation of a simple binary tree of integers mapping to integers
#ifndef _BINTREE
#define _BINTREE

typedef int bintree_key;
typedef int bintree_value;

struct bintree
{
  struct bintree_node* nodes;	//Array of nodes
  struct bintree_node* root;	//The root
};

struct bintree_node
{
  bintree_key key;
  bintree_value value;
  struct bintree_node *left, *right;
};

///This creates a binary tree from a \b sorted! list
void bintree_create_from_ordered_list(struct bintree* tree, int* keys, int* values, int num_entries);

///This destroys the binary tree
void bintree_destroy(struct bintree* tree);

///This gets the element with key \a key
///\return 0 if the element was not found, 1 if it was found
///\note If the element was not found, \a value will not be written to
int bintree_get(struct bintree* tree, bintree_key key, bintree_value* value);

#endif