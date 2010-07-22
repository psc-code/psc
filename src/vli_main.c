
#include <mpi.h>

#include "psc.h"
#include "util/params.h"

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  psc_create("fortran", "fortran", "fortran", "fortran");

  psc_init(NULL);
  psc_integrate();

  MPI_Finalize();
}
