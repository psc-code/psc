
#include <vector_private.h>

#include <stdlib.h>
#include <assert.h>

struct vector_double {
  double *elements;
};

// ----------------------------------------------------------------------
// vector_double_setup

static void
vector_double_setup(struct vector *vec)
{
  struct vector_double *dbl = mrc_to_subobj(vec, struct vector_double);

  dbl->elements = calloc(vec->nr_elements, sizeof(*dbl->elements));
}

// ----------------------------------------------------------------------
// vector_double_destroy

static void
vector_double_destroy(struct vector *vec)
{
  struct vector_double *dbl = mrc_to_subobj(vec, struct vector_double);

  free(dbl->elements);
  dbl->elements = NULL; // just to be safe
}

// ----------------------------------------------------------------------
// vector_double_set_element

static void
vector_double_set_element(struct vector *vec, int i, double val)
{
  struct vector_double *dbl = mrc_to_subobj(vec, struct vector_double);

  assert(dbl->elements);
  assert(i >= 0 && i < vec->nr_elements);

  dbl->elements[i] = val;
}

// ----------------------------------------------------------------------
// vector_double_get_element

static double
vector_double_get_element(struct vector *vec, int i)
{
  struct vector_double *dbl = mrc_to_subobj(vec, struct vector_double);

  assert(dbl->elements);
  assert(i >= 0 && i < vec->nr_elements);

  return dbl->elements[i];
}

// ----------------------------------------------------------------------

struct vector_ops vector_double_ops = {
  .name                  = "double",
  .size                  = sizeof(struct vector_double),
  .setup                 = vector_double_setup,
  .destroy               = vector_double_destroy,
  .set_element           = vector_double_set_element,
  .get_element           = vector_double_get_element,
};
