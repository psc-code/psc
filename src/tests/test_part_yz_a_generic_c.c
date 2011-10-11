
#include "psc_testing.h"
#include "psc_push_particles.h"

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

void
psc_testing_push_yz_a(struct psc *psc, const char *s_push_particles)
{
  if (opt_checks_verbose) {
    mprintf("=== testing push_part() %s\n", s_push_particles);
  }

  psc_testing_dump(psc, s_push_particles);
  psc_push_particles_push_yz_a(psc->push_particles, psc->particles, psc->flds);
  psc_testing_dump(psc, s_push_particles);
}

// ----------------------------------------------------------------------
// check push_particles "cuda_1st" against C "1st" ref

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  struct psc *psc = psc_testing_create_test_yz("fortran", 0, "c");
  psc_setup(psc);
  psc_testing_push_yz_a(psc, "fortran");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz("generic_c", 0, "c");
  psc_setup(psc);
  psc_testing_push_yz_a(psc, "generic_c");
  psc_check_particles_ref(psc, psc->particles, 1e-7, "generic_c");
  psc_destroy(psc);

  psc_testing_finalize();
}
