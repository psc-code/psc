
#include "psc_randomize_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_randomize_fortran_run

static void
psc_randomize_fortran_run(struct psc_randomize *randomize,
			  mparticles_base_t *particles_base)
{
  assert(ppsc->nr_patches == 1);
  
  mparticles_fortran_t *particles = psc_mparticles_get_fortran(particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("fort_randomize", 1., 0, 0);
  }
  prof_start(pr);
  PIC_randomize(psc_mparticles_get_patch_fortran(particles, 0));
  prof_stop(pr);

  psc_mparticles_put_fortran(particles, particles_base);
}

// ======================================================================
// psc_randomize: subclass "fortran"

struct psc_randomize_ops psc_randomize_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_randomize_fortran_run,
};
