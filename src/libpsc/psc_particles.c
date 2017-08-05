
#include "psc.h"
#include "psc_particles_private.h"

#include <mrc_profile.h>
#include <string.h>

// ======================================================================
// psc_particles_init

static void
psc_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_single_by_block_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_cuda_ops);
#endif
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_cuda2_ops);
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_acc_ops);
#endif
}

// ======================================================================
// psc_particles class

struct mrc_class_psc_particles mrc_class_psc_particles = {
  .name             = "psc_particles",
  .size             = sizeof(struct psc_particles),
  .init             = psc_particles_init,
};

