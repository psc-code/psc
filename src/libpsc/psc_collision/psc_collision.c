
#include "psc_collision_private.h"

// ======================================================================
// forward to subclass

void
psc_collision_run(struct psc_collision *collision, mparticles_base_t *particles)
{
  struct psc_collision_ops *ops = psc_collision_ops(collision);
  assert(ops->run);
  ops->run(collision, particles);
}

// ======================================================================
// psc_collision_init

static void
psc_collision_init()
{
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_fortran_ops);
}

// ======================================================================
// psc_collision class

struct mrc_class_psc_collision mrc_class_psc_collision = {
  .name             = "psc_collision",
  .size             = sizeof(struct psc_collision),
  .init             = psc_collision_init,
};

