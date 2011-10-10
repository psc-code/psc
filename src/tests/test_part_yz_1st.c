
#include "psc_testing.h"

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

// ----------------------------------------------------------------------
// check push_particles C "1st" against "generic_c" (C 2nd order) ref
// 
// since the fields are linear, 1st vs 2nd order should not make a
// difference and we expect the same answer

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  struct psc *psc = psc_testing_create_test_yz("1st", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "1st");
  psc_check_continuity(psc, psc->particles, psc->flds, 1e-12);
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz("generic_c", 0, "c");
  psc_setup(psc);
  psc_testing_push_particles(psc, "generic_c");
  psc_testing_push_particles_check(psc, 1e-7, 1e-0);
  psc_destroy(psc);

  psc_testing_finalize();
}
