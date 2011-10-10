
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_sort.h"
#include "psc_randomize.h"
#include "psc_moments.h"
#include "psc_bnd.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

static unsigned int mask;

static bool check_currents = true;
static bool check_particles = true;

// ======================================================================
// psc_test_ops

struct psc_ops psc_test_ops = {
  .name             = "test",
  .size             = sizeof(struct psc_test),
  .create           = psc_test_create,
  .init_field       = psc_test_init_field_linear,
  .init_npt         = psc_test_init_npt_rest,
  .step             = psc_test_step,
};

static struct psc *
create_test(const char *s_push_particles)
{
  struct psc *psc = psc_create(MPI_COMM_WORLD);
  psc->domain.gdims[0] = 1; // make yz
  psc_push_particles_set_type(psc->push_particles, s_push_particles);
  psc_sort_set_type(psc->sort, "countsort2");
  psc_randomize_set_type(psc->randomize, "c");
  psc_moments_set_type(psc->moments, "1st");
  psc_set_from_options(psc);
  psc_setup(psc);
  psc_sort_set_param_int(ppsc->sort, "mask", mask);
  psc_randomize_run(psc->randomize, psc->particles);
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  psc_sort_run(psc->sort, psc->particles);

  return psc;
}

static void
psc_testing_save_ref(struct psc *psc)
{
  psc_save_particles_ref(psc, psc->particles);
  psc_save_fields_ref(psc, psc->flds);
}

static void
psc_testing_check(struct psc *psc, double eps_particles, double eps_fields)
{
  psc_check_continuity(psc, psc->particles, psc->flds, eps_fields);
  if (check_particles) {
    psc_check_particles_ref(psc, psc->particles, eps_particles, "push_part_yz()");
  }
  if (check_currents) {
    psc_check_currents_ref(psc, psc->flds, eps_fields);
  }
}

static void
run_test(struct psc *psc, const char *s_push_particles)
{
  printf("=== testing push_part_yz() %s\n", s_push_particles);

  psc_testing_dump(psc, s_push_particles);
  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  psc_sort_run(psc->sort, psc->particles);
  psc_testing_dump(psc, s_push_particles);
}

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  mrc_params_get_option_bool("check_currents", &check_currents);
  mrc_params_get_option_bool("check_particles", &check_particles);

  // ----------------------------------------------------------------------
  // push_yz 1st order

  struct psc *psc = create_test("1st");
  run_test(psc, "1st");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  // since the fields are linear functions of position, 1st order / 2nd order
  // field interpolation should give the same result
  psc = create_test("generic_c");
  run_test(psc, "generic_c");
  psc_testing_check(psc, 1e-7, 1e-0);
  psc_destroy(psc);

#ifdef xUSE_CUDA
  psc = create_test("cuda_1st");
  run_test(psc, "cuda_1st");
  psc_testing_check(psc, 1e-6, 1e-3);
  psc_destroy(psc);
#endif

  psc_testing_finalize();
}
