
#include "psc.h"
#include "psc_case.h"

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_from_options(_case);
  psc_case_setup(_case);
  psc_case_view(_case);
  psc_integrate(&psc);
  psc_case_destroy(_case);

  MPI_Finalize();
}
