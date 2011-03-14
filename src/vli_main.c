
#include <mpi.h>

#include "psc.h"
#include <mrc_params.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  _psc_case = _psc_case_create(MPI_COMM_WORLD);
  _psc_case_set_from_options(_psc_case);
  _psc_case_setup(_psc_case);
  // FIXME -> _psc_case
  if (psc.prm.from_checkpoint) {
    psc_read_checkpoint();
  }
  _psc_case_view(_psc_case);
  psc_integrate();

  MPI_Finalize();
}
