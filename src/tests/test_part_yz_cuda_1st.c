
#include "psc_testing.h"

// ----------------------------------------------------------------------
// check push_particles "cuda_1st" against C "1st" ref

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops_1);

  struct psc *psc = psc_testing_create_test_yz("1st", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "1st");
  psc_testing_save_ref(psc);
  psc_destroy(psc);

  psc = psc_testing_create_test_yz("cuda_1st", 0, "1st");
  psc_setup(psc);
  psc_testing_push_particles(psc, "cuda_1st");
  psc_testing_push_particles_check(psc, 1e-6, 1e-3);
  psc_destroy(psc);

  psc_testing_finalize();
}
