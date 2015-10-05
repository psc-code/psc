#include <stdlib.h>
#include <assert.h>
#include <bintree.h>

static int numelements = 0;
#ifndef NULL
#define NULL 0
#endif

static void addelement(struct bintree *tree, struct bintree_node **root, bintree_key key, bintree_value value)
{
  //Store the element
  struct bintree_node* node = &tree->nodes[numelements++];
  //Initialize node
  node->key = key;
  node->value = value;
  node->left = NULL;
  node->right = NULL;
  
  //Store node at lowest parent
  struct bintree_node* parent = *root;
  struct bintree_node** target = root;
  
  while(parent != NULL)
  {
    if( key < parent->key )
    {
      if(parent->left != NULL)
      {
	parent = parent->left;
	continue;
      }
      else
      {
	target = &parent->left;
	break;
      }
    }
    else if( key > parent->key )
    {
      if(parent->right != NULL)
      {
	parent = parent->right;
	continue;
      }
      else
      {
	target = &parent->right;
	break;
      }
    }
    else
    {
      assert(0); //You tried to insert a value that already exists
    }
  }
  
  *target = node;
}

static void addfromarray(struct bintree* tree, struct bintree_node** root, int* keys, int* values, int lower, int higher)
{
  if(higher < lower) return;
  
  int median = lower + (higher-lower) / 2;
  
  addelement(tree, root, keys[median], values[median]);
  
  //Left side
  addfromarray(tree, &(*root)->left, keys, values, lower, median-1);
  
  //Right side
  addfromarray(tree, &(*root)->right, keys, values, median+1, higher);
}

void bintree_create_from_ordered_list(struct bintree* tree, int *keys, int *values, int num_entries)
{
  tree->nodes = malloc(sizeof(struct bintree_node) * num_entries);
  tree->root = NULL;

  addfromarray(tree, &tree->root, keys, values, 0, num_entries-1);
  
  numelements = 0;
}

void bintree_destroy(struct bintree* tree)
{
  free(tree->nodes);
  tree->root = NULL;
}

int bintree_get(struct bintree* tree, bintree_key key, bintree_value* value)
{
  struct bintree_node* node = tree->root;
  while( node != NULL )
  {
    if(node->key == key)
    {
      *value = node->value;
      return 1;
    }
    else if(key < node->key)
    {
      node = node->left;
    }
    else
    {
      node = node->right;
    }
  }
  return 0;
}
