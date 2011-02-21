
#include "psc_testing.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  
  struct psc_mod_config conf_fortran = {
    .mod_sort = "fortran",
  };
  psc_create_test_xz(&conf_fortran);
  struct psc_mparticles *particles = &psc.particles;

  psc_randomize();
  psc_sort(particles);
  psc_check_particles_sorted(particles);

  psc_destroy();

  prof_print();
  MPI_Finalize();
}
