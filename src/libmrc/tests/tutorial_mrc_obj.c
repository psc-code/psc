
#include <vector.h>
#include <mrc_params.h>
#include <assert.h>

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

  int nr_elements;
  vector_get_param_int(vec, "nr_elements", &nr_elements);

  for (int i = 0; i < nr_elements; i++) {
    vector_set_element(vec, i, 2. * i);
  }

  for (int i = 0; i < nr_elements; i++) {
    assert(vector_get_element(vec, i) == 2. * i);
  }

  vector_destroy(vec);

  MPI_Finalize();
}
