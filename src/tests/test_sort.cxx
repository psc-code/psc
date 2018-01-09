
#include "psc_testing.h"
#include "psc_randomize.h"
#include "psc_sort.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <mpi.h>

int
main(int argc, char **argv)
{
#if 0
  psc_testing_init(&argc, &argv);
  
  struct psc_case *_case = psc_create_test_xz();
  psc_sort_set_type(ppsc->sort, "fortran");
  psc_case_setup(_case);
  struct psc_mparticles *particles = ppsc->particles;
  psc_randomize_run(ppsc->randomize, particles);
  psc_sort_run(ppsc->sort, particles);
  psc_check_particles_sorted(ppsc, particles);
  psc_case_destroy(_case);

  _case = psc_create_test_xz();
  psc_sort_set_type(ppsc->sort, "countsort");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_randomize_run(ppsc->randomize, particles);
  psc_sort_run(ppsc->sort, particles);
  psc_check_particles_sorted(ppsc, particles);
  psc_case_destroy(_case);

  _case = psc_create_test_xz();
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_randomize_run(ppsc->randomize, particles);
  psc_sort_run(ppsc->sort, particles);
  psc_check_particles_sorted(ppsc, particles);
  psc_case_destroy(_case);

  psc_testing_finalize();
#endif
}
