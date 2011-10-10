
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
// check push_particles "cuda_1vb" against C "1vb" ref

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  struct psc *psc = psc_testing_create_test_yz("1vb", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "1vb");
  psc_check_continuity(psc, psc->particles, psc->flds, 1e-12);
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz("cuda_1vb", 15, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "cuda_1vb");
  psc_testing_push_particles_check(psc, 1e-6, 1e-3);
  psc_destroy(psc);

  psc_testing_finalize();
}
