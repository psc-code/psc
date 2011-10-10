
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

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  // ----------------------------------------------------------------------
  // push_yz 1st order

  struct psc *psc = psc_testing_create_test_yz("1st", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "1st");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  // since the fields are linear functions of position, 1st order / 2nd order
  // field interpolation should give the same result
  psc = psc_testing_create_test_yz("generic_c", 0, "c");
  psc_setup(psc);
  psc_testing_push_particles(psc, "generic_c");
  psc_testing_push_particles_check(psc, 1e-7, 1e-0);
  psc_destroy(psc);

#ifdef xUSE_CUDA
  psc = psc_testing_create_test_yz("cuda_1st", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "cuda_1st");
  psc_testing_push_particles_check(psc, 1e-6, 1e-3);
  psc_destroy(psc);
#endif

  psc_testing_finalize();
}
