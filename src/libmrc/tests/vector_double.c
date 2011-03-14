
#include <vector_private.h>

#include <stdlib.h>
#include <assert.h>

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
// vector_double_set_element

static void
vector_double_set_element(struct vector *vec, int i, double val)
{
  assert(vec->elements);
  assert(i >= 0 && i < vec->nr_elements);

  vec->elements[i] = val;
}

// ----------------------------------------------------------------------
// vector_double_get_element

static double
vector_double_get_element(struct vector *vec, int i)
{
  assert(vec->elements);
  assert(i >= 0 && i < vec->nr_elements);

  return vec->elements[i];
}

// ----------------------------------------------------------------------

struct vector_ops vector_double_ops = {
  .name                  = "double",
  .setup                 = vector_double_setup,
  .destroy               = vector_double_destroy,
  .set_element           = vector_double_set_element,
  .get_element           = vector_double_get_element,
};
