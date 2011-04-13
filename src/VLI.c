
#include "psc.h"
#include "psc_case.h"

#include <mrc_params.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_from_options(_case);
  psc_case_setup(_case);
  psc_case_view(_case);
  psc_case_integrate(_case);
  psc_case_destroy(_case);

  libmrc_params_finalize();
  MPI_Finalize();
}
