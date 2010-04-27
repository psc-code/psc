
#include <mpi.h>

#include "psc.h"
#include "util/params.h"

#define VLA_main_F77 F77_FUNC_(vla_main, VLA_MAIN)

void VLA_main_F77(void);

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);
  
  psc_create("fortran");
  VLA_main_F77();

  MPI_Finalize();
}
