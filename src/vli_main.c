
#include <mpi.h>

#include "psc.h"
#include "util/params.h"

#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)
#define INIT_param_fortran_F77 F77_FUNC_(init_param_fortran, INIT_PARAM_FORTRAN)
#define VLI_main_F77 F77_FUNC_(vli_main, VLI_MAIN)

void INIT_basic_F77(void);
void INIT_param_fortran_F77(void);
void VLI_main_F77(void);

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  psc_create("fortran", "fortran", "fortran");

  psc_init_param();

  SET_param_domain();
  SET_param_psc();
  SET_param_coeff();
  INIT_basic_F77();
  INIT_param_fortran_F77();

  int n_part;
  psc_init_partition(&n_part);
  SET_subdomain();
  psc.f_part = ALLOC_particles(n_part);

  VLI_main_F77();

  MPI_Finalize();
}
