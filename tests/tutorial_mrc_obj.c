
#include <mrc_obj.h>
#include <mrc_params.h>
#include <assert.h>

// ======================================================================
// vector class

struct vector {
  struct mrc_obj obj;
  int nr_elements;
};

MRC_CLASS_DECLARE(vector, struct vector);

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
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct vector *vec = vector_create(MPI_COMM_WORLD);

  assert(vec->nr_elements == 1);
  vec->nr_elements = 3;
  vector_view(vec);

  vector_destroy(vec);

  MPI_Finalize();
}
