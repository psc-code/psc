
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
  struct psc_case *_case = psc_create_test_xz(&conf_fortran);
  psc_case_setup(_case);
  mparticles_base_t *particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  struct psc_mod_config conf_countsort = {
    .mod_sort = "countsort",
  };
  _case = psc_create_test_xz(&conf_countsort);
  psc_case_setup(_case);
  particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  struct psc_mod_config conf_countsort2 = {
    .mod_sort = "countsort2",
  };
  _case = psc_create_test_xz(&conf_countsort2);
  psc_case_setup(_case);
  particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  prof_print();
  MPI_Finalize();
}
