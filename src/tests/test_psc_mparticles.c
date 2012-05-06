
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
  struct psc *psc = psc_testing_create_test_yz("fortran", 0, "1st");
  psc_setup(psc);
  psc_testing_save_ref(psc);
  struct psc_mparticles *mprts;
  if (strcmp(s_type, "c") == 0) {
    mprts = psc_mparticles_get_c(psc->particles, flags);
    psc_mparticles_put_c(mprts, psc->particles, 0);
  } else if (strcmp(s_type, "fortran") == 0) {
    mprts = psc_mparticles_get_fortran(psc->particles, flags);
    psc_mparticles_put_fortran(mprts, psc->particles, 0);
#ifdef USE_CUDA
  } else if (strcmp(s_type, "cuda") == 0) {
    mprts = psc_mparticles_get_cuda(psc->particles, flags);
    psc_mparticles_put_cuda(mprts, psc->particles, 0);
#endif
  } else {
    assert(0);
  }
  psc_check_particles_ref(psc, psc->particles, eps_particles, "mparticles");
  psc_destroy(psc);

  psc_testing_finalize();
}
