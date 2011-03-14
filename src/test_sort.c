
#include "psc_testing.h"
#include "psc_randomize.h"
#include "psc_sort.h"

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
  mparticles_base_t *particles = &psc.particles;

  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(particles);

  psc_destroy(&psc);

  prof_print();
  MPI_Finalize();
}
