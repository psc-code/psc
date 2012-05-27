
#include "psc_randomize_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_randomize_fortran_run

static void
psc_randomize_fortran_run(struct psc_randomize *randomize,
			  struct psc_particles *prts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_randomize", 1., 0, 0);
  }

  struct psc_particles *prts = psc_particles_get_as(prts_base, "fortran", 0);

  prof_start(pr);
  PIC_randomize(prts);
  prof_stop(pr);

  psc_particles_put_as(prts, prts_base, 0);
}

// ======================================================================
// psc_randomize: subclass "fortran"

struct psc_randomize_ops psc_randomize_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_randomize_fortran_run,
};
