
#include <vector_private.h>

#include <stdlib.h>
#include <assert.h>

struct vector_float {
  float *elements;
};

// ----------------------------------------------------------------------
// vector_float_setup

static void
vector_float_setup(struct vector *vec)
{
  struct vector_float *dbl = mrc_to_subobj(vec, struct vector_float);

  dbl->elements = calloc(vec->nr_elements, sizeof(*dbl->elements));
}

// ----------------------------------------------------------------------
// vector_float_destroy

static void
vector_float_destroy(struct vector *vec)
{
  struct vector_float *dbl = mrc_to_subobj(vec, struct vector_float);

  free(dbl->elements);
  dbl->elements = NULL; // just to be safe
}

// ----------------------------------------------------------------------
// vector_float_set_element

static void
vector_float_set_element(struct vector *vec, int i, double val)
{
  struct vector_float *dbl = mrc_to_subobj(vec, struct vector_float);

  assert(dbl->elements);
  assert(i >= 0 && i < vec->nr_elements);

  dbl->elements[i] = val;
}

// ----------------------------------------------------------------------
// vector_float_get_element

static double
vector_float_get_element(struct vector *vec, int i)
{
  struct vector_float *dbl = mrc_to_subobj(vec, struct vector_float);

  assert(dbl->elements);
  assert(i >= 0 && i < vec->nr_elements);

  return dbl->elements[i];
}

// ----------------------------------------------------------------------

struct vector_ops vector_float_ops = {
  .name                  = "float",
  .size                  = sizeof(struct vector_float),
  .setup                 = vector_float_setup,
  .destroy               = vector_float_destroy,
  .set_element           = vector_float_set_element,
  .get_element           = vector_float_get_element,
};
