
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include "psc_bnd.h"
#include "psc_moments.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <mpi.h>

static struct psc_case *
create_test(const char *s_push_particles)
{
  struct psc_case *_case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, s_push_particles);
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_sort_set_param_int3(ppsc->sort, "blocksize", (int [3]) { 1, 8, 8 }); // FIXME
  psc_case_setup(_case);
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

// ----------------------------------------------------------------------
// check push_particles_push_yz against "fortran" ref

// psc_push_particles_type to be tested
static const char *s_type = "fortran";
// threshold for particles
static double eps_particles = 1e-7;
// threshold for fields
static double eps_fields = 1e-7;

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_params_get_option_string("type", &s_type);
  mrc_params_get_option_double("eps_particles", &eps_particles);
  mrc_params_get_option_double("eps_particles", &eps_fields);

  struct psc_case *_case = create_test("fortran");
  psc_testing_push_particles(ppsc, "fortran");
  psc_testing_save_ref(ppsc);
  psc_case_destroy(_case);

  _case = create_test(s_type);
  psc_testing_push_particles(ppsc, s_type);
  psc_testing_push_particles_check(ppsc, eps_particles, eps_fields);
  psc_case_destroy(_case);

  psc_testing_finalize();
}
