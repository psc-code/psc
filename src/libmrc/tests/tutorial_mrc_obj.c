
#include <mrc_obj.h>
#include <mrc_params.h>
#include <assert.h>
#include <stdlib.h>

// ======================================================================
// vector class

struct vector {
  struct mrc_obj obj;
  int nr_elements;
  int *elements;
};

MRC_CLASS_DECLARE(vector, struct vector);

// ======================================================================

// ----------------------------------------------------------------------
// vector_setup

void
_vector_setup(struct vector *vec)
{
  vec->elements = calloc(vec->nr_elements, sizeof(*vec->elements));
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
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct vector *vec = vector_create(MPI_COMM_WORLD);
  vector_set_from_options(vec);
  vector_setup(vec);
  vector_view(vec);
  vector_destroy(vec);

  MPI_Finalize();
}
