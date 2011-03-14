
#include <vector_private.h>

#include <mrc_params.h>
#include <assert.h>

// ----------------------------------------------------------------------
// forwarding to subclass

void
vector_set_element(struct vector *vec, int i, double val)
{
  vector_ops(vec)->set_element(vec, i, val);
}

double
vector_get_element(struct vector *vec, int i)
{
  return vector_ops(vec)->get_element(vec, i);
}

// ----------------------------------------------------------------------
// vector_init

static void
vector_init()
{
  mrc_class_register_subclass(&mrc_class_vector, &vector_double_ops);
  mrc_class_register_subclass(&mrc_class_vector, &vector_float_ops);
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
  .init             = vector_init,
  .param_descr      = vector_descr,
};

