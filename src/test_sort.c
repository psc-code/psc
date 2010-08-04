
#include "psc_testing.h"
#include "util/profile.h"
#include "util/params.h"

#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);
  
  struct psc_mod_config conf_fortran = {
    .mod_sort = "fortran",
  };
  psc_create_test_xz(&conf_fortran);

  PIC_find_cell_indices();
  PIC_randomize();
  psc_sort();
  psc_check_particles_sorted();

  psc_destroy();

  prof_print();
  MPI_Finalize();
}
