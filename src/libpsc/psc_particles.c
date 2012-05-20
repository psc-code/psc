
#include "psc.h"
#include "psc_particles_private.h"

// ======================================================================
// psc_particles_init

static void
psc_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_double_ops);
}

// ======================================================================
// psc_particles class

struct mrc_class_psc_particles mrc_class_psc_particles = {
  .name             = "psc_particles",
  .size             = sizeof(struct psc_particles),
  .init             = psc_particles_init,
};

