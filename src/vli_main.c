
#include <mpi.h>

#include "config.h"

#define VLI_main_F77 F77_FUNC_(vli_main, VLI_MAIN)

void VLI_main_F77(void);

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  
  VLI_main_F77();

  MPI_Finalize();
}
