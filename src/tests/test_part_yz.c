
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

static bool check_currents = true;

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

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields)
{
  printf("=== testing push_part_yz() %s %s\n", s_push_particles,
	 is_ref ? "(ref)" : "");

  struct psc_case *_case = create_test(s_push_particles);
  psc_testing_dump(ppsc, s_push_particles);
  psc_push_particles_run(ppsc->push_particles, ppsc->particles, ppsc->flds);

  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  psc_testing_dump(ppsc, s_push_particles);

  psc_check_continuity(ppsc, ppsc->particles, ppsc->flds, 1e-14);
  if (is_ref) {
    psc_save_particles_ref(ppsc, ppsc->particles);
    psc_save_fields_ref(ppsc, ppsc->flds);
  } else {
    psc_check_particles_ref(ppsc, ppsc->particles, eps_particles, "push_part_yz()");
    if (check_currents) { // only check currents for full pusher
      psc_check_currents_ref(ppsc, ppsc->flds, eps_fields);
    }
  }
  psc_case_destroy(_case);
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

  mrc_params_get_option_bool("check_currents", &check_currents);

  mrc_params_get_option_string("type", &s_type);
  mrc_params_get_option_double("eps_particles", &eps_particles);
  mrc_params_get_option_double("eps_particles", &eps_fields);

  run_test(true, "fortran", 0., 0.);
  run_test(false, s_type, eps_particles, eps_fields);

  psc_testing_finalize();
}
