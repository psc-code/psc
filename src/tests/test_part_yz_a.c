
#include "psc_testing.h"
#include "psc_push_particles.h"

#include <mrc_params.h>

void
psc_testing_push_yz_a(struct psc *psc, const char *s_push_particles)
{
  if (opt_checks_verbose) {
    mprintf("=== testing push_part_yz_a() %s\n", s_push_particles);
  }

  psc_testing_dump(psc, s_push_particles);
  psc_push_particles_push_yz_a(psc->push_particles, psc->particles, psc->flds);
  psc_testing_dump(psc, s_push_particles);
}

// ----------------------------------------------------------------------
// check push_particles_push_yz_a against "fortran" ref

// psc_push_particles_type to be tested
static const char *s_type = "fortran";
// threshold for particles
static double eps_particles = 1e-7;

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_params_get_option_string("type", &s_type);
  mrc_params_get_option_double("eps_particles", &eps_particles);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops_1);

  struct psc *psc = psc_testing_create_test_yz("fortran", 0, "c");
  psc_setup(psc);
  psc_testing_push_yz_a(psc, "fortran");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz(s_type, 0, "c");
  psc_setup(psc);
  psc_testing_push_yz_a(psc, s_type);
  psc_check_particles_ref(psc, psc->particles, eps_particles, s_type);
  psc_destroy(psc);

  psc_testing_finalize();
}
