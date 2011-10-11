
#include "psc_testing.h"
#include <mrc_params.h>

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
  mrc_params_get_option_double("eps_fields", &eps_fields);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops_1);

  struct psc *psc = psc_testing_create_test_yz("fortran", 0, "c");
  psc_setup(psc);
  psc_testing_push_particles(psc, "fortran");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz(s_type, 0, "c");
  psc_setup(psc);
  psc_testing_push_particles(psc, s_type);
  psc_testing_push_particles_check(psc, eps_particles, eps_fields);
  psc_destroy(psc);

  psc_testing_finalize();
}
