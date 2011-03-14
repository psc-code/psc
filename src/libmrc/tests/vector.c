
#include <vector.h>

#include <mrc_params.h>
#include <assert.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// vector_setup

void
_vector_setup(struct vector *vec)
{
  vec->elements = calloc(vec->nr_elements, sizeof(*vec->elements));
}

// ----------------------------------------------------------------------
// vector_destroy

void
_vector_destroy(struct vector *vec)
{
  free(vec->elements);
  vec->elements = NULL; // just to be safe
}

// ----------------------------------------------------------------------
// vector_set_element

void
vector_set_element(struct vector *vec, int i, double val)
{
  assert(vec->elements);
  assert(i >= 0 && i < vec->nr_elements);

  vec->elements[i] = val;
}

// ----------------------------------------------------------------------
// vector_get_element

double
vector_get_element(struct vector *vec, int i)
{
  assert(vec->elements);
  assert(i >= 0 && i < vec->nr_elements);

  return vec->elements[i];
}

// ----------------------------------------------------------------------

#define VAR(x) (void *)offsetof(struct vector, x)
static struct param vector_descr[] = {
  { "nr_elements"       , VAR(nr_elements)      , PARAM_INT(1)  },
  {},
};
#undef VAR

struct mrc_class_vector mrc_class_vector = {
  .name             = "vector",
  .size             = sizeof(struct vector),
  .param_descr      = vector_descr,
  .setup            = _vector_setup,
  .destroy          = _vector_destroy,
};

