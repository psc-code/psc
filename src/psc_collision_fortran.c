
#include "psc_collision_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_collision_fortran_run

static void
psc_collision_fortran_run(struct psc_collision *collision,
			  mparticles_base_t *particles_base)
{
  mparticles_fortran_t particles;
  particles_fortran_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_collision", 1., 0, 0);
  }
  prof_start(pr);
  PIC_bin_coll(&particles.p[0]);
  prof_stop(pr);

  particles_fortran_put(&particles, particles_base);
}

// ======================================================================
// psc_collision: subclass "fortran"

struct psc_collision_ops psc_collision_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_collision_fortran_run,
};
