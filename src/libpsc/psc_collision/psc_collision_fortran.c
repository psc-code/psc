
#include "psc_collision_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_collision_fortran_run

static void
psc_collision_fortran_run(struct psc_collision *collision,
			  mparticles_base_t *particles_base)
{
  mparticles_fortran_t *particles = psc_mparticles_get_fortran(particles_base, 0);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_collision", 1., 0, 0);
  }
  prof_start(pr);
  particles_fortran_t *pp = psc_mparticles_get_patch_fortran(particles, 0);
  PIC_bin_coll(pp);
  prof_stop(pr);

  psc_mparticles_put_fortran(particles, particles_base, 0);
}

// ======================================================================
// psc_collision: subclass "fortran"

struct psc_collision_ops psc_collision_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_collision_fortran_run,
};
