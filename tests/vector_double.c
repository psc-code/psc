
#include <vector_private.h>

#include <stdlib.h>

// ----------------------------------------------------------------------
// vector_double_setup

static void
vector_double_setup(struct vector *vec)
{
  vec->elements = calloc(vec->nr_elements, sizeof(*vec->elements));
}

// ----------------------------------------------------------------------
// vector_double_destroy

static void
vector_double_destroy(struct vector *vec)
{
  free(vec->elements);
  vec->elements = NULL; // just to be safe
}

// ----------------------------------------------------------------------

struct vector_ops vector_double_ops = {
  .name                  = "double",
  .setup                 = vector_double_setup,
  .destroy               = vector_double_destroy,
};
