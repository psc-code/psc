
#include <mpi.h>

#include "psc.h"
#include <mrc_params.h>

#define VLA_main_F77 F77_FUNC_(vla_main, VLA_MAIN)

void VLA_main_F77(void);

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  
  struct psc_mod_config conf = {
  };
  psc_create(&conf);
  VLA_main_F77();

  MPI_Finalize();
}
