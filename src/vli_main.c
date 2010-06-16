
#include <mpi.h>

#include "psc.h"
#include "util/params.h"

#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)
#define VLI_main_F77 F77_FUNC_(vli_main, VLI_MAIN)

void INIT_basic_F77(void);
void VLI_main_F77(void);

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  psc_create("fortran", "fortran", "fortran");

  INIT_basic_F77();
  VLI_main_F77();

  MPI_Finalize();
}
