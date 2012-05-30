
#include "psc_collision_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_collision_c_run

static void
psc_collision_c_run(struct psc_collision *collision,
		    struct psc_particles *prts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("collision", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);

  prof_start(pr);
  PIC_bin_coll(prts);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
}

// ======================================================================
// psc_collision: subclass "c"

struct psc_collision_ops psc_collision_c_ops = {
  .name                  = "c",
  .run                   = psc_collision_c_run,
};
