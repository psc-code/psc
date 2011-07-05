
#include "psc.h"
#include "psc_case.h"

#include <mrc_params.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int checkpoint_step = 0, checkpoint_nmax = 0;
  mrc_params_get_option_int("checkpoint_step", &checkpoint_step);
  mrc_params_get_option_int("checkpoint_nmax", &checkpoint_nmax);
  
  struct psc *psc = psc_read_checkpoint(MPI_COMM_WORLD, checkpoint_step);
  ppsc = psc;
  if (checkpoint_nmax) {
    psc->prm.nmax = checkpoint_nmax;
  }
  psc_view(psc);
  psc_integrate(psc);
  psc_destroy(psc);

  libmrc_params_finalize();
  MPI_Finalize();
}
