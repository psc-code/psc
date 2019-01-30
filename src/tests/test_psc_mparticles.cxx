
#include "psc_testing.h"
#include <mrc_params.h>

#include <string.h>

// ----------------------------------------------------------------------
// check push_particles_push_yz against "fortran" ref

// psc_mparticles type to be tested
static const char *s_type = "fortran";
// threshold for particles
static double eps_particles = 1e-7;
// flags
static unsigned int flags = 0;

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_params_get_option_string("type", &s_type);
  mrc_params_get_option_double("eps_particles", &eps_particles);
  mrc_params_get_option_uint("flags", &flags);

  mrc_class_register_subclass(&mrc_class_psc, (int *) &psc_test_ops_1);

  // FIXME moments, at least
  struct psc *psc = psc_testing_create_test_yz("fortran", 0);// moments "1st";
  psc_setup(psc);
  psc_testing_save_ref(psc);

  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, s_type, flags);
  psc_mparticles_put_as(mprts, psc->particles, 0);

  psc_check_particles_ref(psc, psc->particles, eps_particles, "mparticles");
  psc_destroy(psc);

  psc_testing_finalize();
}
