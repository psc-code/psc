
#include "psc.h"
#include "util/profile.h"

#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  psc_create_test_1("fortran");

  for (int i = 0; i < psc.n_part; i++) {
    psc.f_part[i].cni = random();
  }
  psc_sort();
  psc_check_particles_sorted();

  psc_destroy();

  prof_print();
  MPI_Finalize();
}
