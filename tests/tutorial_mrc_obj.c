
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

struct mrc_class_vector mrc_class_vector = {
  .name             = "vector",
  .size             = sizeof(struct vector),
};

// ======================================================================

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct vector *vec = vector_create(MPI_COMM_WORLD);
  vec->nr_elements = 3;
  assert(vec->nr_elements == 3);
  vector_destroy(vec);

  MPI_Finalize();
}
