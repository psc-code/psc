
#include "psc_collision_private.h"

// ======================================================================
// forward to subclass

void
psc_collision_run(struct psc_collision *collision, mparticles_base_t *particles)
{
  static int st_time_collision;
  if (!st_time_collision) {
    st_time_collision = psc_stats_register("time collision");
  }

  psc_stats_start(st_time_collision);
  struct psc_collision_ops *ops = psc_collision_ops(collision);
  assert(ops->run);
  for (int p = 0; p < particles->nr_patches; p++) {
    ops->run(collision, psc_mparticles_get_patch(particles, p));
  }
  psc_stats_stop(st_time_collision);
}

// ======================================================================
// psc_collision_init

static void
psc_collision_init()
{
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_single_ops);
#ifdef USE_FORTRAN
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_fortran_ops);
#endif
}

// ======================================================================
// psc_collision class

struct mrc_class_psc_collision mrc_class_psc_collision = {
  .name             = "psc_collision",
  .size             = sizeof(struct psc_collision),
  .init             = psc_collision_init,
};

