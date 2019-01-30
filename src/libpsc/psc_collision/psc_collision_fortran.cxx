
#include "psc_collision_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_collision_fortran_run

static void
psc_collision_fortran_run(struct psc_collision *collision,
			  struct psc_mparticles *mprts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_collision", 1., 0, 0);
  }

  PscMparticlesFortran mprts = mprts_base->get_as<PscMparticlesFortran>();

  prof_start(pr);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    PIC_bin_coll(prts);
  }

  prof_stop(pr);

  mprts.put_as(mprts_base);
}

// ======================================================================
// psc_collision: subclass "fortran"

struct psc_collision_ops psc_collision_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_collision_fortran_run,
};
