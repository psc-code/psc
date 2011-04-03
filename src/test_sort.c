
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
  
  struct psc_case *_case = psc_create_test_xz();
  psc_sort_set_type(psc.sort, "fortran");
  psc_case_setup(_case);
  mparticles_base_t *particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  _case = psc_create_test_xz();
  psc_sort_set_type(psc.sort, "countsort");
  psc_case_setup(_case);
  particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  _case = psc_create_test_xz();
  psc_sort_set_type(psc.sort, "countsort2");
  psc_case_setup(_case);
  particles = &psc.particles;
  psc_randomize_run(psc.randomize, particles);
  psc_sort_run(psc.sort, particles);
  psc_check_particles_sorted(&psc, particles);
  psc_case_destroy(_case);

  prof_print();
  MPI_Finalize();
}
