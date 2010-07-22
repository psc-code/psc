
#include <mpi.h>

#include "psc.h"
#include "util/params.h"

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  struct psc_mod_config conf = {
  };
  psc_create(&conf);
  psc_init(NULL);
  psc_integrate();

  MPI_Finalize();
}
